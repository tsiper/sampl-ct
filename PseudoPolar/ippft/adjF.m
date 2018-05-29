% function x=adjF(y)
%
% Adjoint of the operator F
%
% Parameters:
%      y        vector of length n+1 (n even)
%      x        vector of length n
%
% Yoel Shkolnisky 30/03/03

function x=adjF(y)

% Check that the vector y is of length n+1 with n even.
n = length(y)-1;
if mod(n,2)~=0
    error('y must be of length n+1 (n even)');
end
m=2*n+1;

y2 = zeros(1,m);
y2(1:2:m) = y;
x = m.*icfft(y2);
x = x(n/2+1:3*n/2);
