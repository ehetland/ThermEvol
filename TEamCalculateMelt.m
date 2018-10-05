function [x,T,X,mat,MATERIAL] = TEamCalculateMelt(t,xo,T,Xo,mat,...
					  MESH,MATERIAL,AdvectPos)
%
%
% AdvectPos = [bottom_nde_pos top_nde_pos re-injection-point]
%
global ofid fcname version Scales


% calculate partial melt
X = MeltingEq(T,MATERIAL.XTM,MATERIAL.XTB,MATERIAL.XTT);
% ensure the surface does not melt
X(end-1:end) = 0;

if max(MATERIAL.meltcomp)>0
  x = zeros(size(xo));
  dX = X-Xo;
  % calculate compaction due to melting/solidification
  N = size(xo,1);
  TotalCompact = 0;
  % use downward propegation of compaction
  x(N) = xo(N);
  for j=N-1:-1:2
    ThisCompact = -0.5*(xo(j+1)-xo(j-1))*...
      MATERIAL.meltcomp(j)*dX(j);
    x(j) = xo(j) + TotalCompact + ThisCompact/2;
    TotalCompact = TotalCompact+ThisCompact;
  end
  x(1) = xo(1) + TotalCompact;
else
  x = xo;
end

if isempty(AdvectPos)
  return
end

if MATERIAL.Advect.Method==10
  %
  % simply removes a prescribed temperature at a prescribed time
  %
  if t>MATERIAL.Advect.Time & MATERIAL.Advect.Pass==1
    MATERIAL.Advect.Pass = 0;
    apos = [AdvectPos(1):AdvectPos(2)];
    %apos = find(x>=AdvectPos(1)&x<=AdvectPos(2));
    T(apos) = T(apos) - MATERIAL.Advect.Temp./Scales.Temperature;;
    fprintf(ofid,['%s: removing %8.2fC from sill at ',...
	       'time %0.2f yrs in current cycle\n'],...
	    fcname,MATERIAL.Advect.Temp,t.*Scales.Time/3.155e7);
  end
elseif MATERIAL.Advect.Method==11
  if t>MATERIAL.Advect.Time & MATERIAL.Advect.Pass==1
    %
    % removes a prescribed percent of melt (calculates the temp
    % drop from the melt enthalpy, regardless of whether it is
    % there), and does not change the melting equilibrium equation
    %
    MATERIAL.Advect.Pass = 0;
    apos = [AdvectPos(1):AdvectPos(2)];
    %apos = find(x>=AdvectPos(1)&x<=AdvectPos(2));
    % enthalpy in J/kg
    Cmelt = MATERIAL.ME(apos,1) .* T(apos) .* ...
	    Scales.Temperature + MATERIAL.ME(apos,2);
    Xf = max(0,X(apos)-MATERIAL.Advect.Percent);
    Tsink = (Cmelt./MATERIAL.Cp(apos)).*...
	    MATERIAL.rhoM(apos).*MATERIAL.Advect.Percent./...
	    (MATERIAL.rhoM(apos).*Xf + (1-X(apos)).* ...
	     MATERIAL.rhoS(apos));
    T(apos) = T(apos) - Tsink./Scales.Temperature;;
    X(apos) = Xf;
    fprintf(ofid,['%s: removing %6.4f percent melt from sill at ',...
	       'time %0.2f yrs in current cycle (ave temp drop = %8.2fC)\n'],...
	    fcname,MATERIAL.Advect.Percent,...
	    t.*Scales.Time/3.155e7,mean(Tsink));
  end
elseif MATERIAL.Advect.Method==12
  %
  % removes the remaining melt in the sill once the melt falls
  % below a certain value, and re-sets the melting equilibrium
  % curve so that the solidus temp is at the current temp, then
  % contracts the layer by the percent of melt removed
  %
  apos = [AdvectPos(1):AdvectPos(2)];
  %apos = find(x>=AdvectPos(1)&x<=AdvectPos(2));
  if length(AdvectPos)<3
    % did not specify reinjection point
    AdvectPos(3) = NaN;
  end
  dx = (x(apos+1)-x(apos-1))./2;
  TotalMeltPercent = sum(X(apos).*dx)/sum(dx);
  if TotalMeltPercent<=MATERIAL.Advect.Percent & MATERIAL.Advect.Pass==1
    MATERIAL.Advect.Pass = 0;
    fprintf(ofid,['%s: ADVECTION  at time %0.2f yrs\n'],fcname,...
	    t.*Scales.Time/3.155e7);
    % enthalpy in J/kg
    Cmelt = MATERIAL.ME(apos,1) .* T(apos) .* ...
	    Scales.Temperature + MATERIAL.ME(apos,2);
    Tsink = (Cmelt./MATERIAL.Cp(apos)).*...
	    (MATERIAL.rhoM(apos).*X(apos))./...
	    ((1-X(apos)).*MATERIAL.rhoS(apos));
    TotalTsink = sum(Tsink.*dx)/sum(dx);
    % set the melt temp to the depth integrated current temperature
    % of the sill
    MeltTemp = sum(T(apos).*dx)/sum(dx);
    Tshift = T(apos) - MATERIAL.XTT(apos,1);
    %figure(1)
    %clf
    %plot(x,T,'k-o')
    %hold on
    T(apos) = T(apos)-Tsink./Scales.Temperature;
    X(apos) = 0;
    %plot(x,T,'b-o')
    fprintf(ofid,['%s: removing the remaining melt from sill at ',...
	       'time %0.2f yrs in current cycle (depth integrated ',...
	       'temp drop = %8.2fC)\n'],fcname,...
	    t.*Scales.Time/3.155e7,TotalTsink);
    fprintf(ofid,['%s: increasing the solidus/etc temps of the' ...
	       ' sill to the present temp\n'],fcname);
    OldSillXTT = MATERIAL.XTT(apos,1:end);
    OldSillXTM = MATERIAL.XTM(apos,1:end);
    OldSillXTB = MATERIAL.XTB(apos,1:end);
    MATERIAL.XTT(apos,1:end) = MATERIAL.XTT(apos,1:end) + ...
	repmat(Tshift,1,size(MATERIAL.XTT,2));
    MATERIAL.XTB(apos,1:end) = MATERIAL.XTB(apos,1:end) - ...
	MATERIAL.XTM(apos,1:end).*repmat(Tshift,1, ...
					 size(MATERIAL.XTM,2));
    % collapse sill by TotalMeltPercent
    MeltColHeight = sum(dx)*TotalMeltPercent;
    fprintf(ofid,'%s: collapsing sill by %7.2f m\n', ...
	    fcname,MeltColHeight.*Scales.Length);
    [x,MATERIAL] = TEamCollapseSill(x,MATERIAL,apos,MeltColHeight);
    %figure(1)
    %plot(x,T,'g-o')
    if ~isnan(AdvectPos(3))
      % inject collapsed thickness
      [tmp,SillPos] = min(abs(x-AdvectPos(3)));
      NewSill = [AdvectPos(3)-MeltColHeight AdvectPos(3)];
      fprintf(ofid,['%s: injecting %7.2f m of ',...
		 'melt at %9.4f km depth and %8.2fC\n'],...
	      fcname,MeltColHeight.*Scales.Length,...
	      AdvectPos(3).*Scales.Length/1e3,...
	      MeltTemp*Scales.Temperature-273.15);
      % for now just inherit the thermal properties of host
      % material
      SillMat.k = MATERIAL.k(SillPos);
      SillMat.L = MATERIAL.L(SillPos);
      SillMat.Cp = MATERIAL.Cp(SillPos);
      SillMat.rhoS = MATERIAL.rhoS(SillPos);
      SillMat.rhoM = MATERIAL.rhoM(SillPos);
      SillMat.ME = MATERIAL.ME(SillPos,1:end);
      % inherit the melt eq curve for the basalt, but change the
      % liquidus temp to the temp of the melt
      Tshift = MeltTemp - mean(OldSillXTT(:,end),1);
      SillMat.XTT = mean(OldSillXTT,1) + ...
	  repmat(Tshift,1,size(MATERIAL.XTT,2));
      SillMat.XTB  = mean(OldSillXTB,1) - ...
	  mean(OldSillXTM,1).*repmat(Tshift,1,size(MATERIAL.XTM,2));
					      
      SillMat.XTM = mean(OldSillXTM,1);
      % keep the melt compaction rate of the basalt (THIS SHOULD BE
      % FIXED TO BE THE COMPACTION RATE OF THE RESIDUAL MATERIAL!!!)
      SillMat.meltcomp = MATERIAL.meltcomp(apos(1),1:end);
      % inject the new sill
      [x,T,X,mat,MATERIAL,InjPos] = TEamInjectSill(x,T,X,MESH, ...
						   MATERIAL,NewSill, ...
						   MeltTemp,SillMat, ...
						   x(AdvectPos(1:2)));
      %figure(1)
      %plot(x,T,'r-x');
    end
  end
end

% THIS IS OLD STUFF
if 1==0
  METHOD = MATERIAL.Advect.Method;
  if sum(MATERIAL.AdvectMax)>0
    % IF METHOD 1, DOES NOT USE Advect.Layer IN INPUT
    if METHOD==1
      % find each element that the melt fraction exceeds value
      pos = find(X>=MATERIAL.AdvectCriticalMelt & ...
		 MATERIAL.AdvectNum<MATERIAL.AdvectMax & ...
		 MATERIAL.Advect>0);
      if ~isempty(pos)
	% enthalpy in J/kg
	Cmelt = MATERIAL.ME(pos,1).*Tn(pos).*...
		Scales.Temperature+ ...
		MATERIAL.ME(pos,2);
	Tsink(pos) = -1.*(MATERIAL.rhoM(pos)./MATERIAL.rhoS(pos)).*...
	    (X(pos)./(1-X(pos))).*...
	    (Cmelt./MATERIAL.Cp(pos));
	% this is to hack in that temp sink is not point value
	for ip=1:length(pos)
	  Tsink(pos(ip)-1) = Tsink(pos(ip));
	  Tsink(pos(ip)+1) = Tsink(pos(ip));
	end
	fprintf(ofid,['%s: advecting melt at time %0.2f yrs ',...
		   'in current cycle\n'],fcname,...
		t.*Scales.Time/3.155e7)
	X(pos) = 0;
	MATERIAL.AdvectNum(pos) = MATERIAL.AdvectNum(pos)+1;
	mps = find(MATERIAL.AdvectNum>=MATERIAL.AdvectMax);
	MATERIAL.Advect(mps) = 0;
      end
    elseif METHOD==2
      keyboard
    elseif METHOD==3
      % find the average melt fraction in the layers that are turned
      % on
      for k=1:max(MATERIAL.Advect)
	% find melt every MATERIAL.Advect.Layer(k)
	pos = find(MATERIAL.Advect==k);
	if max(X(pos))>=max(MATERIAL.AdvectCriticalMelt(pos))
	  Cmelt = MATERIAL.ME(pos,1).*Tn(pos)+ ...
		  MATERIAL.ME(pos,2);
	  Tsink(pos) = (MATERIAL.rhoM(pos)./MATERIAL.rhoS(pos)).*...
	      (X(pos)./(1-X(pos))).*...
	      (Cmelt./MATERIAL.Cp(pos));
	  X(pos) = 0;
	  MATERIAL.AdvectNum(pos) = MATERIAL.AdvectNum(pos)+1;
	  mps = find(MATERIAL.AdvectNum(pos)>=MATERIAL.AdvectMax(pos));
	  MATERIAL.Advect(pos(mps)) = 0;
	end
      end
      
    end
  end
end



return

