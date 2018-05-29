% function [res1,res2]= slowPPFT(im);
%
% Compute the pseudo-polar Fourier transform directly.
% The computation requires O(n^4) operations.
%
% im    The image whose pseudo-polar Fourier transform should be computed.
%       Must be of a dyadic square size (2^k x 2^k).
%
% Returns res1 and res2 (of size 2n+1xn+1) that contain the pseudo-polar Fourier 
% transform of the input image im.
% res1 contains the values of PPI1(k,l). res2 contains the values of PPI2(k,l). 
% The first argument of Res1 and Res2 corresponds to pseudo-radius and the second argument 
% corresponds to pseudo-angle.
%
% See PPFT.m for more documentation.
%
% See thesis' final PDF p.42 for more information.
%
% Yoel Shkolnisky 22/10/01

function [res1,res2] = slowPPFT(im);

im = flipud(im);
s=size(im);
if (length(s)~=2) %not a 2D image
   error('Input must be a 2D image');
end

if (s(1)~=s(2))
   error('Input image must be square');
end

if (mod(s(1),2)~=0)
   error('Input image must have even side');
end

n=s(1);
ppRows = 2*n+1;
ppCols = n+1;
res1 = zeros(ppRows,ppCols);
res2 = zeros(ppRows,ppCols);

% computation of Res1
for k=-n:n
   for l=-n/2:n/2
      res1(toUnaliasedIdx(k,ppRows),toUnaliasedIdx(l,ppCols))=trigPoly(im,-2*l*k/n,k);
   end
end

% computation of Res2
for k=-n:n
   for l=-n/2:n/2
      res2(toUnaliasedIdx(k,ppRows),toUnaliasedIdx(l,ppCols))=trigPoly(im,k,-2*l*k/n);
   end
end


%%%%%%%%%%%%%%%     Subroutines   %%%%%%%%%%%%%%%%%%
function w=trigPoly(im,xi1,xi2)
% The function evaluates the trigonometric polynomial
%                          n/2    n/2
%       TrigPolyI(xi1,xi2) = sum    sum I(u,v)*exp(-2*pi*i*(xi1*u+xi2*v)/m)       m=2n+1
%                         u=-n/2 v=-n/2
% where:
%       xi1,xi2 -   The points in the frequency domain at which we evaluate the trigonometric polynomial
%       im          -  The 2-D signal (image) for which the trigonometric polynomial is evaluated
%
% This function is used as an auxiliary function for computing the pseudo-polar Fourier transform
% directly (according to the definition of the pseudo-polar Fourier transform).
% The calling function (slowPPFT) uses this auxiliary function to compute the trigonometric 
% polynomial on the pseudo-polar grid points.


s=size(im);
n=s(1);
m=2*n+1;
acc = 0;
for u=lowIdx(n):hiIdx(n)     % x direction
   for v=lowIdx(n):hiIdx(n)  % y direction 
      acc = acc + im(toUnaliasedIdx(v,n),toUnaliasedIdx(u,n))*exp(-2*pi*i*(xi1*u+xi2*v)/m);
   end
end
w=acc;


%Revision record
% 9/12/02   Yoel Shkolnisky     Added comments to the function "trigPoly"
