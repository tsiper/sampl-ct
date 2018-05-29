function y = meyeraux(x)
%MEYERAUX auxiliary Meyer function
% Compute Meyer wavelet auxiliary function 
% v(x) = 35*x^4 - 84*x^5 + 70*x^6 - 20*x^7
% at points x.
%
% INPUT:
%  x                (vector) grid points
%
% OUTPUT:
%  y				(vector) values at given points x
%
%--------------------------------------------------------------------------
% Sören Häuser ~ FFST ~ 2014-07-22 ~ last edited: 2014-07-22 (Sören Häuser)

% % Auxiliary function values.
% p = [-20 70 -84 35 0 0 0 0];
% y = polyval(p,x);
% 
% %y = 0 for x<0, y=1 for x>1
% int1 = y.*(x>=0).*(x<=1);
% y = int1 + (x>1);

% combined to one line for performance
y = polyval([-20 70 -84 35 0 0 0 0],x).*(x>=0).*(x<=1) + (x>1);
