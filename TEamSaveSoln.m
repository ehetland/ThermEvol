function [fcnt,ecnt,tcnt] = TEamSaveSoln(...
    to,xo,To,Xo,mato,SOLVER,MATERIAL,pcnt,fcnt,ecnt,tcnt);
					      
global ofid fcname version Scales

fprintf(ofid,['%s: saving model output at time  ',...
	      '%10.2f kyrs\n'],fcname,...
	to(end).*Scales.Time/3.155e10);

% save full solution when needed
if ~isempty(intersect(SOLVER.CyclesOut,pcnt))
  fcnt=fcnt+1;
  t = to;
  x = xo;
  T = To;
  X = Xo;
  mat = mato;
  SaveFile = sprintf('%s/%s-%04d-%04d.mat',...
		     SOLVER.SolnDir,SOLVER.FileBase,pcnt,fcnt);
  fprintf(ofid,['%s: saving full solution to %s\n'],fcname,SaveFile);
  save(SaveFile,'x','t','T','X','mat','MATERIAL');
end
% save integrated melt in CRUST materials
if SOLVER.TrackCrustMelt==1
    
end

% save the tracked nodes
if ~isempty(SOLVER.TrackNodes)
  ecnt=ecnt+1;
  t = to;
  if SOLVER.TrackNodes(1)>0
    x = xo(SOLVER.TrackNodes,:);
    T = To(SOLVER.TrackNodes,:);
    X = Xo(SOLVER.TrackNodes,:);
    mat = mato(SOLVER.TrackNodes,:);
    fprintf(ofid,['%s: saving solution at tracked nodes ', ...
		  'to %s\n'],fcname,SaveFile);
  else
    x = zeros(length(SOLVER.TrackNodes),size(xo,2));
    T = zeros(length(SOLVER.TrackNodes),size(To,2));
    X = zeros(length(SOLVER.TrackNodes),size(Xo,2));
    mat = abs(SOLVER.TrackNodes);
    for ik=1:length(t)
      for lk=1:length(SOLVER.TrackNodes)
	ps = find(mato(:,ik)==abs(SOLVER.TrackNodes(lk)));
	x(lk,ik) = mean(xo(ps,ik));
	T(lk,ik) = mean(To(ps,ik));
	X(lk,ik) = mean(Xo(ps,ik));
      end
    end
  end
  SaveFile = sprintf('%s/%s_track-%04d-%04d.mat',...
		     SOLVER.SolnDir,SOLVER.FileBase,pcnt,ecnt);  
  save(SaveFile,'x','t','T','X','mat');
end
% save solution at specific times
if ~isempty(SOLVER.OutputSteps)
  tcnt=tcnt+1;
  % ensure last time step is saved
  OutputSteps = unique([SOLVER.OutputSteps(find(...
      SOLVER.OutputSteps<=length(to))),length(to)]);
  t = to(OutputSteps);
  x = xo(:,OutputSteps);
  T = To(:,OutputSteps);
  X = Xo(:,OutputSteps);
  mat = mato(:,OutputSteps);
  SaveFile = sprintf('%s/%s_time-%04d-%04d.mat',...
		     SOLVER.SolnDir,SOLVER.FileBase,pcnt,tcnt);
  fprintf(ofid,['%s: saving solution at specified steps ',...
		'to %s\n'],fcname,SaveFile);
  save(SaveFile,'x','t','T','X','mat');
end
  
  
return

