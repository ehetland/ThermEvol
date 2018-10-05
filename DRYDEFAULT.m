clear

% Thickness          50 m
% Rate               10 kyr
% MinDeltaX          5 m
% Duration           2 Myr
% Initial Range      10-15 km
% B.Temperature      1150 C
% Pressure           3 kbar
% Basalt             DRY
% Granitoid          DRY
% Latent Heat        Even


%%%% PARAMTER DEFINITIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = 'DRYDEFAULT';

Num=600;
MinDeltaX=5;
Cycles=100;

FileBase=sprintf('%s_%d_%0.2f',model,Num,MinDeltaX);
% DIRECTORY NAME BELOW FOR OUTPUT %
SolnDir='/cmld/data4/calogero/TestThermEvol/RunOutput/DRYDEFAULT';

SillThick=50;
Period=10;

Moho=-35.0;
Granitoid=-32.0;
BDT=-15.0;

% SCALING
Scales.Length=1e3;
Scales.Time=1e3*3.155e7;
Scales.Temperature=1e3;

% BOUNDARIES [base moho base top surface]
MAT.Border=[-45 Moho Granitoid BDT 0].*1e3;
% ADVECTION IS CURRENTLY: off
MAT.Advect.Max = zeros(1,length(MAT.Border)-1);


%%%% DEFINE SILL EMPLACEMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IF USING AN ALREADY DETERMINED SET OF SILLS
SILL.x=zeros(Cycles,3);
SILL.x(:,2)=SillPositionSave(1:100,2)*1000;

IF GENERATING A NEW SET OF SILLS
for ii=1:Cycles
   SILL.x(ii,2)=unifrnd((-15000)-SillThick*(ii-1), ...
       (-10000),1,1);
end

% fix sill nodes with node grid
SILL.x(:,2)=round(SILL.x(:,2)./MinDeltaX).*MinDeltaX;
SILL.x(:,1)=SILL.x(:,2)-SillThick;
SILL.x(:,3)=0; % no advection


%%%% MATERIAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load MaterialLibrary.mat;

% MELT FRACTION VS TEMPERATURE (0 --> 100)
% double final one if there aren't enough for vertcat
GabbroXTslope = [0.0006 0.0027 0.0023 0.0023];
GabbroXTinter = [-0.6297 -3.0325 -2.4736 -2.4736];
BasaltXTslope = [0.0013 0.0041 0.0023 0.0030];
BasaltXTinter = [-1.2909 -4.2714 -2.2563 -3.2979]; 
UCrustXTslope = [0.0013 0.0022 0.0034 0.0022];
UCrustXTinter = [-1.2808 -2.2534 -3.7021 -2.1676];
MantleXTslope = 1.0e-3.*[1 1 1 1];
MantleXTinter = -1e3.*[1 1 1 1];
LCrustXTslope = 0.0012.*[1 1 1 1];
LCrustXTinter = -1.23.*[1 1 1 1];

% Surrounding Material Properties
MAT.Density = [MantleRhoS  LCrustRhoS UCrustRhoM UCrustRhoM];
MAT.MeltDensity = [MantleRhoM LCrustRhoM 2545 2545];
MAT.SpecificHeat = [MantleCp LCrustCp UCrustCp UCrustCp];
MAT.LatentHeat = [MantleL LCrustL UCrustL UCrustL];
MAT.ThermalConductivity = [Mantlek LCrustk UCrustk UCrustk];
MAT.XT.slope = [MantleXTslope;LCrustXTslope;UCrustXTslope;UCrustXTslope];
MAT.XT.inter = [MantleXTinter;LCrustXTinter;UCrustXTinter;UCrustXTinter];
MAT.Compaction = [MantleMeltCompaction;LCMeltCompaction;...
		  UCMeltCompaction;UCMeltCompaction];
MAT.Advect.Out = zeros(1,length(MAT.Border)-1);
MAT.Advect.Melt = zeros(length(MAT.Border)-1,2);
MAT.MeltEnergy = zeros(length(MAT.Border)-1,2);
MAT.Advect.Method=0;

% Melt Energy in units of J/g (intercept) & J/kg (slope) 1e3
MAT.MeltEnergy(3,1:2)=[1.5228 -774.21].*1e3;

% SILL Material Properties
SILL.Density = BasaltRhoS;
SILL.MeltDensity = BasaltRhoM;
SILL.SpecificHeat = BasaltCp;
SILL.LatentHeat = BasaltL;
SILL.ThermalConductivity = Basaltk;
SILL.MeltEnergy = [0 1088e3];
SILL.slope = BasaltXTslope;
SILL.inter = BasaltXTinter;
SILL.slopextl = GabbroXTslope;
SILL.interxtl = GabbroXTinter;
SILL.T = 1150+273.15;
SILL.Compaction = BasaltMeltCompaction;


%%%% SET GEOTHERM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define upper and lower boundary temperatures (INITIAL GEOTHERM)
Tbndry = polyval(polyfit([-60e3 0],[1200 0]+273.15,1),MAT.Border([1 end]));


%%%% BUILD STRUCTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (can play around with 'gamma')
MESH = struct('Scales',Scales,...
	      'N',Num,...
	      'delta',MinDeltaX,...
	      'BufferSize',5e4,...
	      'smoothing',20,...
	      'gamma',6,...
	      'iter',1);

% (can play around with tolerance)
SOLVER = struct('Tolerance',1e-5,...
		'FileBase',FileBase,...
		'SolnDir',SolnDir,...
		'Cycles',Cycles,...
		'Period',Period,...
		'Monitor',Period/10,...
		'CyclesOut',[25 50 72 100]);
    
    % 'CyclesOut',[1:10 25:5:Cycles]

ThermEvolAM(SOLVER,MESH,MAT,Tbndry,SILL);







