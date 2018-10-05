function [x,T,X,mat,MATERIAL,InjPos] = TestTEamInjectSill(xo,To,Xo,MESH,MATERIAL, ...
					   Sill,SillTemp,SillMat,varargin)
%
% THIS VERSION OF THE CODE IS EXPLICITLY FOR LEEMAN ET AL 2008 STRATEGY
% EVERY SILL LOADED INTO THE SAME POSITION. (MAC 5.24.16)
%
%
%
%
% Num = number of nodes to use
% To,Xo = temp and melt fraction
% Sill = [xmin xmax] = [xbot xtop] : position of injected sill (Sill.x(p,1:2))
% SillTemp : temperature of injected sill (Sill.T(p))
% SillMat = structure of 
%
global ofid fcname version Scales

%-------------------------------------------------------------------
%INSERTED HERE TEMPORARILY TO MAKE EVERYTHING WORK; FROM THERMEVOLAM
Scales = MESH.Scales;
fcname = 'ThermEvolAM';
version = 1.8;
ofid = 1;
%-------------------------------------------------------------------


% the final mesh, temp field, and melt fraction
x = zeros(MESH.N,1);
T = zeros(MESH.N,1);
X = zeros(MESH.N,1);
mat = zeros(size(x));

SillThickness = diff(Sill(1:2));
SillCenter = mean(Sill(1:2));
% determine which layer sill is placed in
SillLayer = find(MATERIAL.Border<mean(Sill(1:2)),1,'last');

% the new computational border will be thicker by the sill
% thickness
% bord = [MATERIAL.Border(1)-SillThickness MATERIAL.Border(end)];
% UPDATED BORD BELOW
% advect all material boundaries below sill middle down
pbnd = find(MATERIAL.Border<SillCenter);
MATERIAL.Border(pbnd) = MATERIAL.Border(pbnd)-SillThickness;

ubnd = find(MATERIAL.Border>Sill(2));
% ADD IN ADDITIONAL BOUNDARY LAYERS FOR THE SILL
% & updated bord
bord=[MATERIAL.Border(pbnd) SillCenter MATERIAL.Border(ubnd)];
MATERIAL.Border=[MATERIAL.Border(pbnd) Sill(1:2) MATERIAL.Border(ubnd)];


% UPDATE MATERIAL PROPERTIES FOR SILL BOUNDARIES
MATERIAL.Density=[MATERIAL.Density(pbnd) SillMat.Density ...
    MATERIAL.Density(pbnd(end)+1:end)];
MATERIAL.MeltDensity=[MATERIAL.MeltDensity(pbnd) SillMat.MeltDensity ...
    MATERIAL.MeltDensity(pbnd(end)+1:end)];
MATERIAL.SpecificHeat=[MATERIAL.SpecificHeat(pbnd) ...
    SillMat.SpecificHeat ...
    MATERIAL.SpecificHeat(pbnd(end)+1:end)];
MATERIAL.LatentHeat=[MATERIAL.LatentHeat(pbnd) ...
    SillMat.LatentHeat ...
    MATERIAL.LatentHeat(pbnd(end)+1:end)];
MATERIAL.ThermalConductivity=[MATERIAL.ThermalConductivity(pbnd) ...
    SillMat.ThermalConductivity ...
    MATERIAL.ThermalConductivity(pbnd(end)+1:end)];
MATERIAL.Compaction=[MATERIAL.Compaction(pbnd); ...
    SillMat.Compaction; MATERIAL.Compaction(pbnd(end)+1:end)];
MATERIAL.MeltEnergy=[MATERIAL.MeltEnergy(pbnd,:); ...
    SillMat.MeltEnergy; MATERIAL.MeltEnergy(pbnd(end)+1:end,:)];
MATERIAL.XT.slope=[MATERIAL.XT.slope(pbnd,:); ...
    SillMat.slope; MATERIAL.XT.slope(pbnd(end)+1:end,:)];
MATERIAL.XT.inter=[MATERIAL.XT.inter(pbnd,:); ...
    SillMat.inter; MATERIAL.XT.inter(pbnd(end)+1:end,:)];
MATERIAL.XT.To=[MATERIAL.XT.To(pbnd,:); ...
    SillMat.To; MATERIAL.XT.To(pbnd(end)+1:end,:)];
% Advection Properties
MATERIAL.Advect.Max=[MATERIAL.Advect.Max(pbnd) 0 ...
    MATERIAL.Advect.Max(pbnd(end)+1:end)];
MATERIAL.Advect.Out=[MATERIAL.Advect.Out(pbnd) 0 ...
    MATERIAL.Advect.Out(pbnd(end)+1:end)];
MATERIAL.Advect.Melt=[MATERIAL.Advect.Melt(pbnd,:); ...
    [0 0]; MATERIAL.Advect.Melt(pbnd(end)+1:end,:)];

N = MESH.N;
M = find(xo>=SillCenter,1,'first');
% split nodes into those above the middle of the sill and those below
% make new computational domain composed of all nodes below/above
% the sill middle shifted down/up by SillThickness/2 and the add
% n nodes into the domain split
if 1==1
  n = ceil(3.*SillThickness/MESH.delta);
  while 1
    x1 = [xo(1:M-1)-1.5*SillThickness;...
	  linspace(SillCenter-1.5*SillThickness,...
		   SillCenter+1.5*SillThickness,n)';...
	  xo(M:N)+1.5*SillThickness];
    % node position at top of Sill on new mesh
    [tmp,M1] = min(abs(x1-Sill(2)));
    % node position at bottom of Sill on new mesh
    [tmp,Q1] = min(abs(x1-Sill(1)));
    SillPos = round((Q1+M1)/2);
    if SillPos-n>N/40
      break
    else
      fprintf(ofid,['%s: WARNING ran out of nodes below the sill, ',...
		    'position decreasing resolution of the mesh ',...
		    'density function in the sill to keep some ', ...
		    'nodes below the sill\n'],fcname);
      n = round(n/2);
      if n<2
	fprintf(ofid,['%s: ERROR could not reconcile sill ',...
		      'position and current mesh, aborting, try' ...
		       ' re-running with more nodes\n'],fcname);
	save DebugOut_TEamInjectSill.mat
	return
      end
    end
  end
else
  % ensure that there is a node at the SillCenter
  n=1;
  x1 = [xo(1:M-1)-0.5*SillThickness;...
	SillCenter;...
	xo(M:N)+0.5*SillThickness];
end
% advect node positions below/above the sill center by SillThickness/2
xop = [xo(1:M-1)-SillThickness/2;xo(M:N)+SillThickness/2];
% map the existing temperature field onto the split nodal
% positions, soothly interpolating the temps over the added sill
% Thickness - this is used for the mesh density dependent on the
% pre-sill temperature field, we use melt fraction field to get
% mesh density in new sill region - use a smooth interpolation so
% that the extrapolation is wel behaved
T1 = interp1(xop,To,x1,'spline','extrap');
% map the melt-fraction onto the new nodal positions
X1 = interp1(xop,Xo,x1,'spline','extrap');
% to ensure density of nodal positions in the neighborhood of the
% sill, set melt fraction to 1 over a 3*SillThickness region
% centered on the sill center
X1(M+1:M+n-2) = 1.0;

if length(varargin)
  % supplied additional hold intervals to not smooth Txx over (this
  % is for the most recent sill
  Hold = varargin{1}';
  [tmp,M2] = min(abs(x1-Hold(2)));
  [tmp,Q2] = min(abs(x1-Hold(1)));
else
  Hold = [];
  M2 = [];
  Q2 = [];
end

% determine mesh density function 
% do not use temp for sill for the gradient to avoid singularities,
% the melt fraction will get the sill location, 
phi = TEamMeshDensity(x1,T1,X1,MESH,sort([Q1-2 M1+2 Q2-2 M2+2]));

AdaptBorders = [bord(1) SillCenter bord(end)];
% AdaptBorders=bord;
AdaptPos = [1 SillPos-n N];
if ~isempty(Hold)
  HoldCenter = mean(Hold);
  HoldPos = round((M2+Q2)/2);
  if (HoldPos-n)>0
    AdaptBorders = [bord(1) HoldCenter SillCenter bord(end)];
    AdaptPos = [1 HoldPos-n SillPos-n N];
    [AdaptPos,spos] = sort(AdaptPos);
    AdaptBorder = AdaptBorders(spos);
  else
    fprintf(ofid,['%S: WARNING ran out of nodes below hold ',...
		  'position, proceeding remeshing without ',...
		  'holding mesh at former sill location\n'],...
	    fcname);
  end
end
% do not map temp/melt field here, not registered on x1, sill
% accomodation done below
save('AdaptMeshInput.mat','x1','MESH','phi','AdaptBorders',...
     'AdaptPos','n');


x = TEamAdaptMesh(x1,[],[],MESH,MATERIAL,phi,AdaptBorders,AdaptPos,n);
if ~isempty(find(isnan(x)))
  save DebugOut_TEamInjectSill.mat
end


% map original temp/melt field onto new mesh & then add sill, this
% ensures that we do not smooth sill during interpolation
Top = interp1(xo,To,x,'linear','extrap');
Xop = interp1(xo,Xo,x,'linear','extrap');
clear To Xo
To = Top;
Xo = Xop;
% node position at top of Sill
[tmp,M] = min(abs(x-Sill(2)));
% node position at bottom of Sill
[tmp,Q] = min(abs(x-Sill(1)));
Q = Q+1;

%PlotDomainTmp
%keyboard

% advect positions down, i.e. set xo to be the positions that To and
% Xo are registered at AFTER sill is injected
xo = x;
xo(1:M) = x(1:M) - SillThickness;

% copy over the temp/melt field above the sill, material properties
% are already set on these nodes
T(M+1:N) = To(M+1:N);
X(M+1:N) = Xo(M+1:N);

% map the advected temp/melt fields below the sill top to the
% node positions below the sill bottom
T(1:Q-1) = interp1(xo(1:M),To(1:M),x(1:Q-1),'linear','extrap');
X(1:Q-1) = interp1(xo(1:M),Xo(1:M),x(1:Q-1),'linear','extrap');

% add in basalt with uniform temperature and pure melt
InjPos = [Q M];
T(Q:M) = SillTemp;
X(Q:M) = 1.0;

fprintf(ofid,['%s: injecting a %8.2fC sill from %9.4f ',...
	      'to %9.4f km in material layer %d, ',...
	      'approximately %7.2f (%7.2f) m thick\n'],...
	fcname,SillTemp.*Scales.Temperature-273.15,...
	x(Q).*Scales.Length/1e3,...
	x(M).*Scales.Length/1e3,SillLayer,...
	diff(x([Q M])).*Scales.Length,...
	diff(x([Q-1 M+1])).*Scales.Length);
	%(mean(x(M:M+1))-mean(x(Q-1:Q))).*Scales.Length);

MATERIAL.k(:) = 0;

% adding material properties of the sill nodes
% props = fieldnames(MATERIAL);
% mat(Q:M) = 0;
% for k=1:length(props)
%   field = getfield(MATERIAL,props{k});
%   if ~isstruct(field) & size(field,1)==N
%     if isempty(strfind(props{k},'Advect'))
%       field(Q:M,:) = repmat(getfield(SillMat,props{k}),M-Q+1,1);
%     else
%       field(Q:M,:) = 0;
%     end
%     MATERIAL = setfield(MATERIAL,props{k},field);
%   end
% end


% set material properties for all non sill nodes on new mesh
acnt=0;
for k=1:length(MATERIAL.Border)-1
   if k==SillLayer
    ps = setdiff(find(MATERIAL.Border(k)<x&x<=MATERIAL.Border(k+ ...
						  1)),[Q:M]);
  else
    ps = find(MATERIAL.Border(k)<x&x<=MATERIAL.Border(k+1));
  end
  if size(ps,1)==1
    ps = ps';
  end
  if k==1
    if isempty(ps)
      keyboard % THIS HAS COME UP, NOT SURE WHY
    end
    ps = unique([1;ps]);
  end
  if (k==length(MATERIAL.Border)-1) & (max(ps)<MESH.N)
   ps = unique([ps;MESH.N]);
  end
  mat(ps) = k;
  MATERIAL.k(ps) = MATERIAL.ThermalConductivity(k);
  MATERIAL.L(ps) = MATERIAL.LatentHeat(k);
  MATERIAL.Cp(ps) = MATERIAL.SpecificHeat(k);
  MATERIAL.rhoS(ps) = MATERIAL.Density(k);
  MATERIAL.rhoM(ps) = MATERIAL.MeltDensity(k);
  % properties used for advection of melt
  if MATERIAL.Advect.Method<=10&MATERIAL.Advect.Method>0
    MATERIAL.AdvectCriticalMelt(ps) = MATERIAL.Advect.Melt(k,1);
    MATERIAL.AdvectResidualMelt(ps) = MATERIAL.Advect.Melt(k,2);
  end
  MATERIAL.ME(ps,1:2) = repmat(MATERIAL.MeltEnergy(k,1:2), ...
			       length(ps),1);
  MATERIAL.meltcomp(ps,1) = MATERIAL.Compaction(k);
  % NOT SURE WHAT TO DO WITH THIS STUFF
  %if MATERIAL.Advect.Out(k)==0
  %  MATERIAL.AdvectMax(ps) = 0;
  %else
  %  BorderTemp = [MATERIAL.Border(k):...
  %		  MATERIAL.Advect.Layer(k)./Scales.Length:...
  %		  MATERIAL.Border(k+1)];
  %  for kk=1:length(BorderTemp)-1
  %    acnt = acnt+1;
  %    pstmp = find(BorderTemp(kk)<x&x<=BorderTemp(kk+1));
  %    MATERIAL.AdvectFlag(pstmp) = acnt;
  %  end
  %  if ~isfield(MATERIAL.Advect,'Max')
  %    MATERIAL.AdvectMax(ps) = 1e10;
  %    fprintf(ofid,['%s: no maximum number of advections for' ...
  %	 ' layer %d specified, using 1e10\n'],...
  %      fcname,k);
  % else
  %   MATERIAL.AdvectMax(ps) = MATERIAL.Advect.Max(k);
  % end
  %end
  % properties used for calculating the melt content at a given temperature
  MATERIAL.XTM(ps,:) = repmat(MATERIAL.XT.slope(k,:),length(ps),1);
  MATERIAL.XTB(ps,:) = repmat(MATERIAL.XT.inter(k,:),length(ps),1);
  MATERIAL.XTT(ps,:) = repmat(MATERIAL.XT.To(k,:),length(ps),1);
end

%MATERIAL.AdvectCriticalMelt = round(MATERIAL.AdvectCriticalMelt.*1e3)./1e3;
%MATERIAL.Advect = round(MATERIAL.Advect);

return
