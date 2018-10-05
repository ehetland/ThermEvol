function ThermEvolAM(SOLVER,MESH,MATERIAL, ...
		   Tbndry,Sill,varargin)

global ofid fcname version Scales

Scales = MESH.Scales;
fcname = 'ThermEvolAM';
version = 1.8;

if isfield(SOLVER,'LogFile')
  ofid = fopen(sprintf('%s.log',SOLVER.LogFile),'w');
else
  ofid = 1;
end

fprintf(ofid,['Starting %s V%0.1f (1D thermal evolution w/ adaptive' ...
	      ' meshing)\n'], fcname,version);

% check the solver structure to make sure everything is there, if
% not set it to defaut value
if ~isfield(SOLVER,'FileBase')
  SOLVER.FileBase = 'ThermEvolAMrun';
end
if ~isfield(SOLVER,'SolnDir')
  SOLVER.SolnDir = 'ThermEvol';
end
if ~isfield(SOLVER,'Tolerance');
  SOLVER.Tolerance = 1e-4;
  fprintf(ofid,'%s: using default solver tolerance of %0.2e\n',...
	  fcname,SOLVER.Tolerance);
end
if ~isfield(SOLVER,'TrackNodes')
  SOLVER.TrackNodes = [];
end
if ~isfield(SOLVER,'TrackCrustMelt')
  SOLVER.TrackCrustMelt = 0;
end
if ~isfield(SOLVER,'OutputSteps')
  SOLVER.OutputSteps = [];
end
if ~isfield(SOLVER,'CyclesOut')
  SOLVER.CyclesOut = [];
end
if ~isfield(SOLVER,'Cycles');
  SOLVER.Cycles = 1;
  fprintf(ofid,['%s: SOLVER.Cycles not set, computing over one', ...
		' period\n'], fcname);
end
if ~isfield(SOLVER,'Period');
  SOLVER.Period = 1;
  fprintf(ofid,['%s: SOLVER.Period not set, computing for', ...
		' period=1\n'], fcname);
end
if ~isfield(SOLVER,'Monitor')
  SOLVER.Monitor = sum(SOLVER.Period)*10*SOLVER.Cycles;
end
if ~isfield(SOLVER,'TmpThreshold')
  SOLVER.TmpThreshold = [];
end
if ~isfield(SOLVER,'TmpTest')
  SOLVER.TmpTest=[];
  if length(SOLVER.TmpThreshold)~=length(SOLVER.TmpTest)
    fprintf(ofid,['%s: Number of Temperature Thresholds does not ', ...
		  ' match number of Temperature Test values\n'], ...
	    fcname);
  end
end
% if ~isfield(SOLVER,'MaxTemp')
%     SOLVER.MaxTemp=[];
% end
% if ~isfield(SOLVER,'MeanTemp')
%     SOLVER.MeanTemp=[];
% end

% check the mesh structure to make sure everything is there, if
% not set it to defaut valu
if ~isfield(MESH,'Scales')
  MESH.Scales = struct('Length',1,'Time',1,'Temperature',1);
  fprintf(ofid,['%s: WARNING characteristic scales are not ',...
		'set in MESH.Scales structure, using ',...
		'characteristic scales of 1 for all fields\n'],...
	  fcname);
end
if ~isfield(MESH,'N');
  MESH.N = 250;
  fprintf(ofid,['%s: WARNING number of total nodes not set in MESH.N, ',...
		'using MESH.N = %d\n',],fcname,MESH.N);
end
if ~isfield(MESH,'delta')
  MESH.delta = 0.1;
  fprintf(ofid,['%s: WARNING minimum nodal spacing not set in' ...
		' MESH.delta, using MESH.delta = %0.2f\n'], ...
	  fcname,MESH.delta);
end
if ~isfield(MESH,'iter')
  MESH.iter = 1;
end
if ~isfield(MESH,'gamma');
  MESH.gamma = 5;
end
if ~isfield(MESH,'smoothing');
  MESH.smoothing = 5
end
if ~isfield(MESH,'BufferSize')
  MESH.BufferSize = 1e4;
end
if ~isfield(MESH,'TargetBoundary')
  MESH.TargetBoundary = 2;
end
if ~isfield(MESH,'TargetFraction')
  MESH.TargetFraction = 0.75;
end




if MATERIAL.Advect.Method==10
  if ~isfield(MATERIAL.Advect,'Temp') | ...
	~isfield(MATERIAL.Advect,'Time')
    fprintf(ofid,['%s: ERROR need field for Temp and Time in',...
	       'Advect structure\n'],fcname);
    return
  end
end
if MATERIAL.Advect.Method==11
  if ~isfield(MATERIAL.Advect,'Percent') | ...
	~isfield(MATERIAL.Advect,'Time')
    fprintf(ofid,['%s: ERROR need field for Percent and Time in',...
	       'Advect structure\n'],fcname);
    return
  end
end
if MATERIAL.Advect.Method==12
  if ~isfield(MATERIAL.Advect,'Percent') 
    fprintf(ofid,['%s: ERROR need field for Percent in',...
	       'Advect structure\n'],fcname);
    return
  end
end

% expand Period to a vector of length Cycles
if length(SOLVER.Period)>1 & length(SOLVER.Period)~=SOLVER.Cycles
  fprintf(ofid,['%s: number of specified cycles does not match ', ...
	     'the given number of periods, using the specified ', ...
	     'periods\n'],fcname);
  SOLVER.Cycles = length(SOLVER.Period);
end
if length(SOLVER.Period)==1 & SOLVER.Cycles>1
  fprintf(ofid,'%s: only one period specified for multiple cycles\n',...
	  fcname);
  SOLVER.Period = repmat(SOLVER.Period,1,SOLVER.Cycles);
end
if size(SOLVER.Period,2)==1
  SOLVER.Period = SOLVER.Period';
end

Tc = [700:1:2000];
for k=1:length(MATERIAL.Border)-1
  [Xc(:,k),MATERIAL.XT.slope(k,:),MATERIAL.XT.inter(k,:), ...
   MATERIAL.XT.To(k,:)] = MeltingEq(Tc,MATERIAL.XT.slope(k,:), ...
				    MATERIAL.XT.inter(k,:));
end
[Xc(:,k+1),Sill.slope,Sill.inter,Sill.To] = ...
      MeltingEq(Tc,Sill.slope,Sill.inter);
  
  % MC 5/4/17
[Xc(:,k+2),Sill.slopextl,Sill.interxtl,Sill.Toxtl] = ...
    MeltingEq(Tc,Sill.slopextl,Sill.interxtl);


% non-dimemnsionalize terms that need to be used non-dimensionally
% below
MESH.delta = MESH.delta/Scales.Length;
MATERIAL.Border = MATERIAL.Border./Scales.Length;
if isfield(Sill,'H')
  Sill.RelativePosition = 1;
  fprintf(ofid,'%s: using relative placement in layer\n',fcname);
  if size(Sill.H,1)<SOLVER.Cycles
    fprintf(ofid,['%s: only one, or not enough, sill thicknesses given for',...
		  ' %d cycles, using the same thickness for all injections\n'],...
	    fcname,SOLVER.Cycles);
    Sill.H = repmat(Sill.H,SOLVER.Cycles,1);
  end
  Sill.H = Sill.H./Scales.Length;
else
  Sill.RelativePosition = 0;
  Sill.x = Sill.x./Scales.Length;
end
Sill.T = Sill.T./Scales.Temperature;
Sill.To = Sill.To./Scales.Temperature;
Sill.Toxtl = Sill.Toxtl./Scales.Temperature; % MC 5/4
Sill.slope = Sill.slope.*Scales.Temperature;
Sill.slopextl = Sill.slopextl.*Scales.Temperature; % MC 5/4
MATERIAL.XT.slope = MATERIAL.XT.slope.*Scales.Temperature;
MATERIAL.XT.To = MATERIAL.XT.To./Scales.Temperature;
Tbndry = Tbndry./Scales.Temperature;

% expand Sill elements to the cycle length
if size(Sill.x,1)<SOLVER.Cycles
  fprintf(ofid,['%s: only one, or not enough, sill location given for',...
	   ' %d cycles, using the same location for all injections\n'],...
	  fcname,SOLVER.Cycles);
  Sill.x = repmat(Sill.x(1,1:end),SOLVER.Cycles,1);
end
if length(Sill.T)<SOLVER.Cycles
  fprintf(ofid,['%s: only one, or not enough, sill temperatures given for',...
	   ' %d cycles, using the same temperature for all injections\n'],...
	  fcname,SOLVER.Cycles);
  Sill.T = repmat(Sill.T(1),SOLVER.Cycles,1);
end



MATERIAL.k = zeros([MESH.N 1]);
MATERIAL.L = zeros([MESH.N 1]);
MATERIAL.Cp = zeros([MESH.N 1]);
MATERIAL.rhoS = zeros([MESH.N 1]);
MATERIAL.rhoM = zeros([MESH.N 1]);
MATERIAL.meltcomp = zeros([MESH.N 1]);
% set each distinct layer as a unique integer flag, so
% MATERIAL.Advect.Layer is only needed along with MATERIAL.Advect.Melt
MATERIAL.AdvectNum = zeros([MESH.N 1]);
MATERIAL.AdvectMax = zeros([MESH.N 1]);
MATERIAL.AdvectCriticalMelt = zeros([MESH.N 1]);
MATERIAL.AdvectResidualMelt = zeros([MESH.N 1]);
MATERIAL.ME = zeros(MESH.N,size(MATERIAL.MeltEnergy,2));
% slopes (XT.M), intercepts (XT.B) and temperatures (XT.T) of the XT
% melting equilibrium curveTemp
MATERIAL.XTM = zeros(MESH.N,size(MATERIAL.XT.slope,2));
MATERIAL.XTB = zeros(MESH.N,size(MATERIAL.XT.inter,2));
MATERIAL.XTT = zeros(MESH.N,size(MATERIAL.XT.To,2));

% rename fields in Sill structure
Sill.k = Sill.ThermalConductivity;
Sill.L = Sill.LatentHeat;
Sill.Cp = Sill.SpecificHeat;
Sill.rhoS = Sill.Density;
Sill.rhoM = Sill.MeltDensity;
Sill.XTM = Sill.slope;
Sill.XTB = Sill.inter;
Sill.XTMxtl = Sill.slopextl; % MC 5/4
Sill.XTBxtl = Sill.interxtl; % MC 5/4
Sill.XTT = Sill.To;
Sill.XTTxtl = Sill.Toxtl; % MC 5/4
Sill.ME = Sill.MeltEnergy;
Sill.meltcomp = Sill.Compaction;


save(sprintf('%s/Cntrl_%s.mat',SOLVER.SolnDir,SOLVER.FileBase), ...
     'SOLVER','MESH','MATERIAL','Tc','Xc', 'Sill','Tbndry')

SillPositionSave=Sill.x;
save(sprintf('%s/Sillx_%s.mat',SOLVER.SolnDir,SOLVER.FileBase), ...
     'SillPositionSave')

if sum(MATERIAL.AdvectResidualMelt)>0
  fprintf(ofid,'%s: residual melt not enabled yet!\n',fcname);
end

TEamCalculate(SOLVER,MESH,MATERIAL,Tbndry,Sill);


return



