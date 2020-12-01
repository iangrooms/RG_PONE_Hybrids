function OUT = RHS(~,Y,p)
% Computes RHS for multiscale L96 model
% Inputs:
%   p   struct with fields K, J, h, and F
%   Y   variables Y_{j,k}
% Output:
%   OUT dY/dt computed from new model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following assumes p.K=41, and N=p.J*p.K
% In retrospect the FFTs would be a lot faster if I didn't use K=41,
% since it's prime.

YHat = fft(Y(:));
XHat = [YHat(1:21);YHat(end-19:end)]*(1/p.J);
X = real(ifft(XHat));

nX = circshift(X,1).*(circshift(X,-1)-circshift(X,2));
NX = interpft(nX,p.J*41);

OUT = p.h*circshift(Y,-1).*(circshift(Y,1)-circshift(Y,-2)) - Y + p.F + NX;