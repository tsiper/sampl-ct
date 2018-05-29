% function im = adjRadon(r1,r2) 
%
% Computes the adjoint 2-D discrete Radon transform.
%
% r1,r2    The sections of the discrete Radon transform whose adjoint should be computed.
%          r1 and r2 must be of size (2n+1)x(n+1) as results from the function Radon.
%
% See also Radon, adjPPFT.
%
% Yoel Shkolnisky 9/2/02

function im = adjRadon(r1,r2)

%Check if the input is of valid size 2x(2n+1)x(n+1)
s1 = size(r1);
s2 = size(r2);

if (sum(s1-s2)~=0)
   error('r1 and r2 must have the size');
end

if (mod(s1(1),2)~=1) | (mod(s1(2),2)~=1)
   error('pp1 and pp2 must be of size (2n+1)x(n+1)');
end

n = [(s1(1)-1)/2;s1(2)-1];
if (n(1)~=n(2))
   error('Input parameter must be of size (2n+1)x(n+1)');
end

n=s1(1); % n - number of rows

pp1 = cfftd(r1,1)/n;
pp2 = cfftd(r2,1)/n;
im = optimizedAdjPPFT(pp1,pp2);

%Revision record
% 15/1/03	Yoel Shkolnisky		Used cfftd instead of column-wise cfft