function [varargout] = MeltingEq(T,M,B,varargin)
%
% assumes a piece-wise linear T-X line between the solidus and
% liquidus
%
% X = MeltingEq(T,M,B,[To],[Xo])
% [X,M,B,To] = MeltingEq(T,M,B,[To],[Xo])
%
% INPUT:
% for N temperatures and p linear segments between TS and TL
% T = temperature [1 N]
% M = slopes of the melting curves [p N]
% B = intercepts of the melting curves [p N]
% OPTIONAL INPUT:
% To = [TS Ti TL], where Ti are the temperatures at the intermediary
% segment breaks (if To is not specified, calculates it)
% Xo = [Xi], where Xi are the melt fractions at the intermediary
% segment breaks (only used if M=B=[])
%
% if B=[] and To is specified, calculates the intercepts from the
% slopes and temperatures at the break-points
%
% if M=B=[] and To and Xo are specfied, calculates the slopes and
% intercepts
%
% if M=B=[] and only To is specfied, calculates a single segment
% based on liquidus to solidus
%
% OUTPUT:
% X = melt fraction at the temperatures T
% OPTIONAL OUTPUT:
% M = slopes of the segments
% B = intercepts of the segments
%

if length(varargin)==0 & ~isempty(B) & ~isempty(M)
  % find the temperatures at the breaks in the linear segments
  To = zeros(size(M,1),size(M,2)+1);
  % solidus
  To(:,1) = -B(:,1)./M(:,1);
  % liquidus
  To(:,end) = (1-B(:,end))./M(:,end);
  % intermediary break points
  for k=2:size(To,2)-1
    To(:,k) = -diff(B(:,k-1:k),1,2)./diff(M(:,k-1:k),1,2);
    To(find(To(:,k)==0),k) = To(find(To(:,k)==0),k-1);
    To(find(isnan(To(:,k))),k) = To(find(isnan(To(:,k))),k-1);
    To(find(isinf(To(:,k))),k) = To(find(isinf(To(:,k))),1);
  end
elseif length(varargin)>=1 & isempty(B) & ~isempty(M)
  % calculate the intercepts
  To = varargin{1};
  B = zeros(size(M));
  B(:,1) = -To(:,1).*M(:,1);
  B(:,end) = 1-To(:,end).*M(:,end);
  for k=2:size(M,2)-1
    B(:,k) = M(:,k-1).*To(:,k) + B(:,k-1) - M(:,k).*To(:,k);
  end
elseif length(varargin)>=1 & isempty(B) & isempty(M)
  % calculate the slopes and intercepts
  To = varargin{1};
  if length(varargin)>=2
  Xo = [zeros(size(varargin{2},1),1) ...
	varargin{2} ...
	ones(size(varargin{2},1),1)];
  else
    Xo = repmat([0 1],size(To,1),1);
  end
  M = zeros(1,length(To)-1);
  B = zeros(size(M));
  for k=1:size(To,2)-1
    M(:,k) = diff(Xo(:,k:k+1),1,2)./diff(To(:,k:k+1),1,2);
    B(:,k) = Xo(:,k)-M(:,k).*To(:,k);
  end
elseif length(varargin)>=1 & ~isempty(B)
  To = varargin{1};
else
  X = zeros(size(T));
  fprintf(1,'MeltingEq: input is not correct\n');
  return
end

X = zeros(size(T));

for k=1:size(M,2)
  X = X + (M(:,k).*T+B(:,k)).*...
      (To(:,k)<=T&T<To(:,k+1));
end
X(find(To(:,k+1)<=T)) = 1;

varargout{1} = X;
if nargout>=3
  varargout{2} = M;
end
if nargout>=3
  varargout{3} = B;
end
if nargout>=4
  varargout{4} = To;
end

return