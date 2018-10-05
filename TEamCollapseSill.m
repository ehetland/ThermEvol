function [x,MATERIAL] = TEamCollapseSill(xo,MATERIAL,spos,Collapse)
%
% (x,spos,Percent)
%
% x = node positions
% spos = positions of nodes to collapse
% Collapse = amount to remove from sill (in column height)
%
global ofid fcname version Scales

x = zeros(size(xo));
N = size(xo,1);

% all positions at the top of sill and above stay the same position
x(spos(end):N) = xo(spos(end):N);
% squeeze positions in sill by collapse height
m = (xo(spos(end))-xo(spos(1))-Collapse)/...
    (xo(spos(end))-xo(spos(1)));
b = xo(spos(end))*(1-m);
x(spos) = m*xo(spos)+b;
% all positions at the bottom of sill and below move up by collapse
% column height
x(1:spos(1)) = xo(1:spos(1))+Collapse;

% advect all material boundaries below sill bottom up
pbnd = find(MATERIAL.Border<min(x(spos)));
MATERIAL.Border(pbnd) = MATERIAL.Border(pbnd)+Collapse;


return


