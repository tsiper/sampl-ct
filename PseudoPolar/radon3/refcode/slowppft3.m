% function pp=slowppft3(im)
%
% Compute the 3-D pseudo-polar Fourier transform directly.
% The computation requires O(n^6) operations.
%
% See ppft3.m for documentation of input and output parameters.
%
% Yoel Shkolnisky 30/01/03

function pp=slowppft3(im)

verifyImage(im);

s=size(im);
n=s(1);
pp=zeros(3,3*n+1,n+1,n+1);

for k=-3*n/2:3*n/2
    for l=-n/2:n/2
        for j=-n/2:n/2
            coord = toUnaliasedCoord([k,l,j],[3*n+1,n+1,n+1]);
            pp(1,coord{:}) = trigPoly(im,k,-2*l*k/n,-2*j*k/n);
            pp(2,coord{:}) = trigPoly(im,-2*l*k/n,k,-2*j*k/n);
            pp(3,coord{:}) = trigPoly(im,-2*l*k/n,-2*j*k/n,k);
        end
    end
end
            

%%%%%%%%%%%%%%%     Subroutines   %%%%%%%%%%%%%%%%%%
function s=trigPoly(im,ox,oy,oz)
% The function evaluates the trigonometric polynomial
%                           n/2-1   n/2-1   n/2-1
%             FI(ox,oy,oz)=  sum     sum     sum  I(u,v,w)exp(-2*pi*i(u*ox+v*oy+w*oz)/m)   m=3n+1
%                           u=-n/2  v=-n/2  w=-n/2
% where:
%       ox,oy,oz   - The points in the frequency domain at which we evaluate the trigonometric polynomial
%       im         - The 3-D to transform to the frequency domain
%
% This function is used as an auxiliary function for computing the pseudo-polar Fourier transform
% directly (according to the definition of the pseudo-polar Fourier transform).
% The calling function (slowppft3) uses this auxiliary function to compute the trigonometric 
% polynomial on the pseudo-polar grid points.

s=size(im);
n=s(1);
m=3*n+1;
acc = 0;
for u=lowIdx(n):hiIdx(n)     % x direction
   for v=lowIdx(n):hiIdx(n)  % y direction 
       for w=lowIdx(n):hiIdx(n)
           coord = toUnaliasedCoord([u,v,w],[n,n,n]);
           acc = acc + im(coord{:})*exp(-2*pi*1i*(u*ox+v*oy+w*oz)/m);
       end
   end
end
s=acc;
