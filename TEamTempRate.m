function dT = TEamTempRate(x,T,X,MATERIAL)

global ofid fcname version Scales


% T(:,j-1),K,B,XT,Tsink

% T(:,j-1) -> T
% T(:,j)   -> Tn
% X(:,j-1) -> X
% X(:,j)   -> Xn
% 

dT = zeros(size(T));
%keyboard
% % 03/28/18 MC Edit to test distribution of Latent Heat.
% for j=1:length(X)
%     if MATERIAL.XTM(j,4)==75.50
%         if X(j)>=0.25
%             MATERIAL.L(j)=0.33*MATERIAL.L(j);
%         elseif X(j)<0.25 && X(j)>0
%             MATERIAL.L(j)=3*MATERIAL.L(j);
%         end
%     end
% end



% use a simple density mixing relation that is for X in volume
% percent melt content
rho = (MATERIAL.rhoM-MATERIAL.rhoS).*X + MATERIAL.rhoS;

% determine the slope values at each point (re-dimensionalize slope
% by temp, as L is dimensional)
mt = zeros(size(dT));
for k=1:size(MATERIAL.XTM,2)
  mt = mt + MATERIAL.XTM(:,k).*(MATERIAL.XTT(:,k)<=T&T<MATERIAL.XTT(:,k+1))/...
       Scales.Temperature;
end


% determine multiplier
Kp = rho.*(MATERIAL.Cp + MATERIAL.L.*mt);
A = (1./Kp).*Scales.Time./Scales.Length.^2;

% find average thermal conductivity at the mid-points
k = MATERIAL.k;
kd = k(1:end-1)+diff(k)./2;

i = [2:length(dT)-1];
dT(i,1) = A(i) .* ...
	  ( ...
	      kd(2:end  ).*(T(i+1)-T(i  ))./(x(i+1)-x(i  )) - ...
	      kd(1:end-1).*(T(i  )-T(i-1))./(x(i  )-x(i-1)) )./ ...  
      (0.5.*(x(i+1)-x(i-1)));


return

