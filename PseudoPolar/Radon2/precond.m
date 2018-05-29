% function Y=precond(X);
%
% Explicit form of the preconditioner used to invert the pseudo-polar Fourier transform
% using the CG method.
% This function is not used directly and should be used as a reference.
% The application of the preconditioner to the adjoint pseudo-polar Fourier transform
% is incorporated into the function PrecondAdjPPFT
%
%   X   A matrix of size (2n+1)x(n+1). The function applies the preconditioner
%       to the matrix X.
%
% Yoel Shkolnisky 17/12/02


function Y=precond(X);

[rows,cols] = size(X);

if (mod(rows,2)~=1) | (mod(cols,2)~=1)
   error('X must be of size (2n+1)x(n+1)');
end

s = [(rows-1)/2;cols-1];
if (s(1)~=s(2))
   error('X must be of size (2n+1)x(n+1)');
end

m = rows; %at this point m=2n+1
n = (m-1)/2;
alpha = 2*(n+1)/(n*m);

pre=zeros(rows,cols);

for k=-n:n
   if k==0
      pre(toUnaliasedIdx(k,rows),:) = 1/(m^2);
   else
      pre(toUnaliasedIdx(k,rows),:) = abs(k*alpha);
   end
end

Y=pre.*X/(rows^2);

%Example:
% a = magic(4);
% [p1,p2] = ppft(a);
% a1 = precond(p1);
% a2 = precond(p2);
% b = adjppft(a1,a2);
% c = ptp(magic(4));
% max(max(abs(b-c)))
% ans =
%  4.4964e-015
