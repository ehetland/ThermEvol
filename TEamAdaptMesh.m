function [varargout] = TEamAdaptMesh(xo,To,Xo,MESH,MATERIAL,phi, ...
    bord,pos,n)

global ofid fcname version Scales

%-------------------------------------------------------------------
%INSERTED HERE TEMPORARILY TO MAKE EVERYTHING WORK; FROM
%THERMEVOLAM
Scales = MESH.Scales;
fcname = 'ThermEvolAM';
version = 1.8;
ofid = 1;
%-------------------------------------------------------------------


x = zeros(pos(end),1);

if isempty(phi) 
  if ~isempty(To) & ~isempty(Xo)
    phi = TEamMeshDensity(xo,To,Xo,MESH);
  else
    fprintf(ofid,['%s: ERROR need to supply either temp & melt ',...
		  'OR mesh density function - not adapting mesh\n'],...
	    fcname);
    varargout{1} = xo;
    return
  end
end

for k=1:length(pos)-1
  % number of nodal differences in segment
  num = pos(k+1)-pos(k);
  A = [1,zeros(1,num);...
       [-eye(num),zeros(num,1)]+[zeros(num,1),eye(num)];...
       zeros(1,num),1];

  L = bord(k+1)-bord(k);
  psi = (L - num*MESH.delta)/sum(phi([pos(k)+1:pos(k+1)]+n));  
  PHI = [bord(k);...
	 phi([pos(k)+1:pos(k+1)]+n).*psi+MESH.delta;...
	 bord(k+1)];
  x(pos(k):pos(k+1)) = inv(A'*A)*A'*PHI;
end

if ~isempty(find(isnan(x)))
  save DebugOut_TEamAdaptMesh.mat
end


% MAKE SURE ALL MATERIAL REGIONS HAVE AT LEAST ONE NODE
% (only works for one sillcenter, needs to be updated later on for
% presence of a secondary sill, not worrying about that right now 5/10)
%x = TEamBorderShift(x,MATERIAL,bord(2:end-1));
% ---------------------------------------------------------------------


fprintf(ofid,['%s: mesh adapted to conditions with min dX=%7.2f m ',...
	      'and max(dx)/min(dx) = %0.4f, ',...
	      'new computation domain is (%9.4f,%9.4f) km\n'],fcname,...
	min(diff(x(2:end-1))).*Scales.Length,...
	max(diff(x))/min(diff(x)),...
	x(1).*Scales.Length./1e3,...
	x(end).*Scales.Length./1e3);

if ~isempty(To)
  fprintf(ofid,['%s: mapping temperature field onto new node positions\n'],...
	  fcname);
  T = interp1(xo,To,x,'linear','extrap');
end
if ~isempty(To)
  fprintf(ofid,['%s: mapping melt fraction field onto new node positions\n'],...
	  fcname);
  X = interp1(xo,Xo,x,'linear','extrap');
end

varargout{1} = x;
if nargout>=2
  varargout{2} = T;
end
if nargout>=3
  varargout{3} = X;
end

return
