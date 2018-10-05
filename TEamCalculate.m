function TEamCalculate(SOLVER,MESH,MATERIAL,Tbndry,Sill)

global ofid fcname version Scales

% initialize solution arrays
t = zeros(1,MESH.BufferSize);
x = zeros(MESH.N,MESH.BufferSize);
T = zeros(MESH.N,MESH.BufferSize);
X = zeros(MESH.N,MESH.BufferSize);
mat = zeros(size(x));

% for initial node positions, use uniform spacing
x(:,1) = [MATERIAL.Border(1):...
	  diff(MATERIAL.Border([1 end]))/(MESH.N-1):...
	       MATERIAL.Border(end)]';



% initial temp profile is linear geotherm
T(:,1) = polyval(polyfit(x([1 end],1),Tbndry',1),x(:,1));

ptime = zeros([SOLVER.Cycles 1]);
ptime(1) = cputime;
delT = 0.5.*MESH.delta^2;
MonitorTime = SOLVER.Monitor;

% TRACKING MAX TEMPERATURE & MEAN TEMPERATURE (MC)
MaxTempinit=max(T(find(x(:,1)<(-10) & x(:,1)>(-15))));
MaxTemp(1)=MaxTempinit;
MeanTempinit=mean(T(find(x(:,1)<(-10) & x(:,1)>(-15))));
MeanTemp(1)=MeanTempinit;
% DETAIL
MaxTempDinit=max(T(find(x(:,1)<(-10) & x(:,1)>(-15))));
MaxTempD(1)=MaxTempDinit;
MeanTempDinit=mean(T(find(x(:,1)<(-10) & x(:,1)>(-15))));
MeanTempD(1)=MeanTempDinit;
ttotal(1)=0;
rcnt=1;
GMTR=zeros(1,11);


for pcnt = 1:SOLVER.Cycles
  
  fprintf(ofid,['%s: START OF PERIOD %d of %d\n'],fcname,pcnt,SOLVER.Cycles);
  
  j = 1;


 
  fprintf(ofid,['%s: adapting mesh to current temperature field ',...
		'before sill injection to ensure %d of %d nodes ',...
		'are in the upper %d layers\n'],...
	  fcname,round(MESH.TargetFraction*MESH.N),MESH.N,...
	  length(MATERIAL.Border)-MESH.TargetBoundary);
  [x(:,j),T(:,j),X(:,j)] = TEamAdaptMesh(x(:,j),T(:,j),X(:,j),...
					 MESH,MATERIAL,[],...
					 MATERIAL.Border([1 MESH.TargetBoundary end]),...
					 [1 round(MESH.N*(1-MESH.TargetFraction)) MESH.N],0);
  
                 
  %figure(1)
  %clf
  %subplot(2,1,1)
  %PlotDX(x(:,1));
  %subplot(2,1,2);
  %plot(x(:,1),T(:,1));
  %keyboard

  if Sill.RelativePosition==1
    Sillx(1:2) = MATERIAL.Border(Sill.x(pcnt,1))+ ...
	diff(MATERIAL.Border(Sill.x(pcnt,1)+[0:1]))*Sill.x(pcnt,2) + ...
	0.5.*[-1 1].*Sill.H(pcnt);
    if size(Sill.x,2)==4
      Sillx(3) = MATERIAL.Border(Sill.x(pcnt,3))+ ...
	diff(MATERIAL.Border(Sill.x(pcnt,3)+[0:1]))*Sill.x(pcnt,4);
    end
  else
    Sillx = Sill.x(pcnt,1:3);
    % MC 5/4/17 %%%%%%%%%%%%%%%%%%%
    if pcnt > 1
        Sillxold=Sill.x(pcnt-1,1:3);
    else
        Sillxold=0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
  fprintf(ofid,['%s: injecting sill from %9.4f to %9.4f km ',...
	   ' at time %0.2f yrs, cycle %d\n'],fcname,...
	  Sillx(1:2).*Scales.Length./1e3,...
	  t(j).*Scales.Time/3.155e7,pcnt); 

  % INJECT INPUT
%   save('InjectSillInput.mat','x','T','X','MESH','MATERIAL', ...
%        'Sillx','Sill');
%   keyboard
%
% REVISED 0416MAC; 

  [x(:,1),T(:,1),X(:,1),mat(:,1),MATERIAL,InjPos]=TestTEamInjectSill(...
      x(:,1),T(:,1),X(:,1),MESH,MATERIAL,...
      Sillx(1:2),Sillxold,Sill.T(pcnt),Sill);
  % MC added Sillxold 5/4/17 %%%%%%%%%%%%%%%%

% InjPos are the most recent injection point of the sill
  if MATERIAL.Advect.Method>=10
    AdvectPos = Sillx(1:end);
    AdvectPos(1:2) = InjPos;
    MATERIAL.Advect.Pass = 1;
  else 
    AdvectPos = [];
  end

  [x(:,1),T(:,1),X(:,1),mat(:,1),MATERIAL] = TEamCalculateMelt(...
      t(j)-sum(SOLVER.Period(1:pcnt-1)),x(:,1),T(:,1),X(:,1),mat(:,1),...
      MESH,MATERIAL,AdvectPos);

  % a counter for output in each period
  fcnt = 0;
  ecnt = 0;
  tcnt = 0;
  %tic
  totalloop=0;
  fid=1;
  while t(j)<sum(SOLVER.Period(1:pcnt))
    j=j+1;
    
    attempt = 0;
    while 1
      attempt = attempt+1;
      if attempt>120;
	fprintf(ofid,'\n\t solution did not converge at step %d\n\n',j);
    keyboard
      end
      
      % convergence does NOT include potential cooling from melt
      % advection that my happen in the time step to be calculated!
      T1 = T(:,j-1);
      dT1 = TEamTempRate(x(:,j-1),T(:,j-1),X(:,j-1),MATERIAL);
        
      T2 = T1+0.5.*delT.*dT1; % TSINK
      dT2 = TEamTempRate(x(:,j-1),T2,X(:,j-1),MATERIAL);

      T3 = T1+delT.*(-dT1+2.*dT2); % TSINK
      dT3 = TEamTempRate(x(:,j-1),T3,X(:,j-1),MATERIAL);
    
      % second order update
      To2 = T1+(delT./2).*(dT1+dT3); % TSINK
        
      % third order update
      To3 = T1+(delT./6).*(dT1+4*dT2+dT3); % TSINK
    
      % difference gives error associated with T2
      eT = (To2-To3);
    
      % check if error is less than tolerance
      error = norm(eT(isfinite(eT)),'inf');
      if error<SOLVER.Tolerance
	break;
      end
      % step is unacceptable, so modify time step
      if 1==0
	fprintf(ofid,...
		['step %d: rejecting attempt %d, err = %0.4e ',...
		 'with dt = %0.2e\n'],...
	j,attempt,error,delT);
      end
      delT = 0.9*delT*(SOLVER.Tolerance/error)^(1/3);
    end
    totalloop=totalloop+attempt;
    % solutions converged at this time step
    if 1==0
      fprintf(ofid,...
	      ['\tstep %d: accepting attempt %d, err = %0.4e ',...
	       'with dt = %0.2e\n'],...
	      j,attempt,error,delT);
    end
    % update time
    t(j) = t(j-1)+delT;

    % accept third order solution as update
    T(:,j) = To3; 
    % calculate the melt at the present time step/temperature
    [x(:,j),T(:,j),X(:,j),mat(:,j),MATERIAL] = TEamCalculateMelt(...
	t(j)-sum(SOLVER.Period(1:pcnt-1)),x(:,j-1),T(:,j),X(:,j-1),mat(:,j-1),...
	MESH,MATERIAL,AdvectPos);
    
    if t(j)>MonitorTime
      MonitorTime = MonitorTime+SOLVER.Monitor;
      fprintf(ofid,['%s: running w/o issue at time   %11.2f kyrs; ',...
		 'Temp = (%8.2f,%8.2f)C, '],fcname,...
	      (t(j)-sum(SOLVER.Period(1:pcnt-1)))*...
	      Scales.Time/3.155e10,...
	      min(T(:,j))*Scales.Temperature-273.15,...
	      max(T(:,j))*Scales.Temperature-273.15);
      fprintf(ofid,['x = (%9.4f,%9.4f)km; dx = (%8.3f,%8.3f)m; ',...
 		 'dt = %9.4f yrs\n'],...
	      x(1,j)*Scales.Length/1e3,...
	      x(end,j)*Scales.Length/1e3,...
	      min(diff(x(:,j)))*Scales.Length,...
	      max(diff(x(:,j)))*Scales.Length,...
	      delT*Scales.Time/3.155e7);
    end
    
    % grow time step
    if error~=0
      delT = 0.9*delT*(SOLVER.Tolerance/error)^(1/3);
    end
    
    %DETAILED SAVE OF MEAN AND MAX T WITH TIME
    Tnow=T(:,j);
    xnow=x(:,j);
    MaxTempDadd=max(Tnow(find(xnow<(-10) & xnow>(-15-0.05*pcnt))));
    MaxTempD=[MaxTempD MaxTempDadd];
    MeanTempDadd=mean(Tnow(find(xnow<(-10) & xnow>(-15-0.05*pcnt))));
    MeanTempD=[MeanTempD MeanTempDadd];
    
    
    
    if t(j)+delT>sum(SOLVER.Period(1:pcnt))
      % if time step is too long, shorten and compute only the 3rd
      % order update
      delT = abs(sum(SOLVER.Period(1:pcnt))-t(j));
      j=j+1;
      t(j) = t(j-1)+delT;
      
      T1 = T(:,j-1);
      dT1 = TEamTempRate(x(:,j-1),T(:,j-1),X(:,j-1),MATERIAL);
      T2 = T1+0.5.*delT.*dT1; % TSINK
      dT2 = TEamTempRate(x(:,j-1),T2,X(:,j-1),MATERIAL);
      T3 = T1+delT.*(-dT1+2.*dT2); % TSINK
      dT3 = TEamTempRate(x(:,j-1),T3,X(:,j-1),MATERIAL);
      To3 = T1+(delT./6).*(dT1+4*dT2+dT3); % TSINK
      T(:,j) = To3;
      % calculate the melt at the present time step/temperature
      [x(:,j),T(:,j),X(:,j),mat(:,j),MATERIAL] = TEamCalculateMelt(...
	  t(j)-sum(SOLVER.Period(1:pcnt-1)),x(:,j-1),T(:,j),X(:,j-1),mat(:,j-1),...
	  MESH,MATERIAL,AdvectPos);
      fprintf(ofid,['%s: calculated solution through time %0.4f kyrs, ',...
		 'period %d of %d\n'],fcname,...
	      t(j)*Scales.Time/3.155e10,pcnt,SOLVER.Cycles);
      
      break;
    end
    if j>=MESH.BufferSize
      %totalloop
      %toc
      [fcnt,ecnt,tcnt] = TEamSaveSoln(t,x,T,X,mat,SOLVER, ...
					     MATERIAL,pcnt,fcnt,ecnt,tcnt);
      totalloop=0;
      %tic

      % re-initialize solution arrays
      t(1) = t(j);
      t(2:end) = 0;
      x(:,1) = x(:,j);
      x(:,2:end) = 0;
      T(:,1) = T(:,j);
      T(:,2:end) = 0;
      mat(:,1) = mat(:,j);
      mat(:,2:end) = 0;
      % reaply boundary conditions
      T(1,:) = T(1,1);
      T(end,:) = T(end,1);
      X(:,1) = X(:,j);
      X(:,2:end) = 0;
      j=1;
    end
  end
  
  if pcnt==SOLVER.Cycles
    tpos = [1:j];
  else
    tpos = [1:j-1];
  end

  %TEMPERATURE THRESHOLD CHECK (MC 1603)
%   [SOLVER.TmpThreshold]=TEamTmpThreshold(...
%     t(tpos),x(:,tpos),T(:,tpos),X(:,tpos),mat(:,tpos),...
%     SOLVER,MATERIAL,pcnt);
  TCh=T(:,tpos(end));
  xch=x(:,tpos(end));
  MaxTempadd=max(TCh(find(xch<(-10) & xch>(-15-0.05*pcnt))));
  MaxTemp(pcnt+1)=MaxTempadd;
  %keyboard
  MeanTempadd=mean(TCh(find(xch<(-10) & xch>(-15-0.05*pcnt))));
  MeanTemp(pcnt+1)=MeanTempadd;
  ttotal=[ttotal;t(1:tpos(end-1))'];
  
  % FOR FINDING MELT AFTER 100 YEARS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % time amount is easily changable, needs to be hardwired in later
  % if we decide to keep this.
  %
  % thundex: t hundred index
  %
  %thundex=find(t(:)>=(((pcnt-1)*SOLVER.Period(pcnt))+0.1),1,'first');
%   if pcnt==1
%       GMTROLD=zeros(1,11);
%   else
%       GMTROLD=GMTR(pcnt-1,:);
%   end
  % MELT THICKNESS SUBROUTINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[GMT(pcnt,[1:11]),GMTR(pcnt,[1:11]),GMTH(pcnt,[1:11]),...
%    BMT(pcnt,[1:11])]=TEamMeltThickness(...
%      MATERIAL,X(:,tpos),x(:,tpos),t(tpos),pcnt,...
%  
  
%   [GMT(pcnt,[1:11]),GMTR(pcnt,[1:11]),GMTH(pcnt,[1:11]),...
%       GMTF(pcnt,[1:11]),GMTHF(pcnt,[1:11]),GMTWH(pcnt,[1:11]),...
%     BMT(pcnt,[1:11])]=TEamMeltThickness(...
%       MATERIAL,X,x,t,pcnt,SOLVER,tpos,GMTROLD);
  
  


%   if pcnt>=100
%   keyboard
%   end
  
  [fcnt,ecnt,tcnt] = TEamSaveSoln(...
      t(tpos),x(:,tpos), T(:,tpos),X(:,tpos),mat(:,tpos),...
      SOLVER,MATERIAL,pcnt,fcnt,ecnt,tcnt);
  ptime(pcnt+1) = cputime;
  fprintf(ofid,'%s: %s to compute over cycle %d\n',...dbcont
	  fcname,TimeString(diff(ptime(pcnt:pcnt+1))),pcnt);
  
  % Breaks if Temperature Threshold is reached (MC 1603)
%   if SOLVER.TmpThreshold(end)==1
%     break
%   end

  % re-initialize solution arrays
  t(1)=t(j);
  t(2:end) = 0;
  x(:,1) = x(:,j);
  x(:,2:end) = 0;
  T(:,1) = T(:,j);
  T(:,2:end) = 0;
  mat(:,1) = mat(:,j);
  mat(:,2:end) = 0;
  % reapply boundary conditions
  T(1,:) = T(1,1);
  T(end,:) = T(end,1);
  X(:,1) = X(:,j);
  X(:,2:end) = 0;
end

SaveTempArray=sprintf('%s/%s_tempstructures.mat',...
    SOLVER.SolnDir,SOLVER.FileBase);
save(SaveTempArray,'MaxTemp','MeanTemp','MaxTempD','MeanTempD','ttotal');

% SaveMeltThick=sprintf('%s/%s_melthick.mat',...
%     SOLVER.SolnDir,SOLVER.FileBase);
% save(SaveMeltThick,'GMT','GMTR','GMTH','GMTF','GMTHF','GMTWH','BMT');

fprintf(ofid,'%s: %s for total computation over %d cycles\n',...
	fcname,TimeString(ptime(end)-ptime(1)),SOLVER.Cycles);

return



function [timestr] = TimeString(timenum)
%
% string = TimeString(time)
%
% forms a text string of the form "xx units' where xx=the time and
% units = seconds, minutes, hours, days
%
% E Hetland, Jan 20087
global ofid fcname version

if timenum<2*60;
  timestr = 'seconds';
elseif timenum<2*60*60
  timenum = timenum/60;
  timestr = 'minutes';
elseif timenum<1.5*24*60*60
  timenum = timenum/3600;
  timestr = 'hours';
else
  timenum = timenum/86400;
  timestr = 'days';
end

timestr = sprintf('%0.3f %s',timenum,timestr);

return

