function phi = TEamMeshDensity(x1,T1,X1,MESH,varargin)
%
% TEamMeshDensity(x,T,X,MESH,[hold])
% sets the rho in hold(i):hold(i+1) i=1,2,4,.. to max, and does not
% smooth
%
global ofid fcname version Scales

N = size(x1,1);
Txx = TEamLaplacian(x1,T1); 
if isfield(MESH,'MaxCurvature')
  Txx(find(abs(Txx)>MESH.MaxCurvature)) = MESH.MaxCurvature;
end
rho = abs(Txx).^(1/4);
rho = rho  + MESH.gamma.*X1;
rho = 1+rho;
rho0 = rho;
%% ensure that the current injected sill has the maximum resolution
%% don't smooth in the sill/hold regions, smooth outside regions up to
%% sill/hold regions...
j = [2:N-1];
if length(varargin)>0
  hold = varargin{1};
  for k=1:2:length(hold)
    rho(hold(k)-1:hold(k+1)+1) = max(rho);
    j = setdiff(j,[hold(k):hold(k+1)]);
  end
end
  
for c=1:MESH.smoothing
  rho(1) = max([rho(1) rho(1)/2 + rho(2)/2]);
  rho(N) = max([rho(N) rho(N-1)/2 + rho(N)/2]);
  %j = [2:N-1];
  rho(j) = rho(j-1)./4 + rho(j)./2 + rho(j+1)./4;
end
phi = 1./rho;
phi = phi-min(phi);
phi = phi./max(phi);

return