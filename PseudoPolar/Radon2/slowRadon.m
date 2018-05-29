%function [res1,res2] = slowRadon(im);
%
% Compute the pseudo-Radon transform directly.
% The computation requires O(n^3) operations.
%
% im    The image whose discrete Radon transform should be computed.
%       Must be real (no imaginary components) and of a dyadic square size (2^k x 2^k).
%
% Returns res1 and res2 (of size 2n+1xn+1) that contain the discrete Radon
% transform of the input image im.
% res1 contains the radon values of basically horizontal lines. res2 contains the Radon
% values of basically vertical lines. 
% The first argument of res1 and res2 corresponds to pseudo-radius and the second argument 
% corresponds to pseudo-angle.
%
% See Radon for more information.
%
% See thesis' final PDF for more information.
%
% Yoel Shkolnisky 22/10/01

function [res1,res2] = slowRadon(im);

im = flipud(im);

s=size(im);
% validity test of the input
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

% computation of res1
for t=-n:n
   for l=-n/2:n/2
      slope = 2*l/n;
      acc   = 0;
      
      for u=-n/2:n/2-1
         acc = acc + I1(im,n,u,slope*u+t);
      end
      res1(toUnaliasedIdx(t,ppRows),toUnaliasedIdx(l,ppCols))=acc;
   end
end     
    
% computation of res2

for t=-n:n
   for l=-n/2:n/2
      slope = 2*l/n;
      acc   = 0;
      
      for v=-n/2:n/2-1
         acc = acc + I2(im,n,slope*v+t,v);
      end
      res2(toUnaliasedIdx(t,ppRows),toUnaliasedIdx(l,ppCols))=acc;
   end
end     


%%%%%%%%%%%%%%%     Subroutines   %%%%%%%%%%%%%%%%%%

% computation of I1 - trigonometric interpolation of I along the columns (along the y-axis)
function z = I1(im,n,u,y)
m=2*n+1;
acc = 0;

for v=-n/2:n/2-1
   acc = acc + im(toUnaliasedIdx(v,n),toUnaliasedIdx(u,n))*dirichlet(y-v,m);
end
z=acc; 
%disp(num2str(z));


% computation of I2 - trigonometric interpolation of I along the columns (along the y-axis)
function z = I2(im,n,x,v)
m=2*n+1;
acc = 0;

for u=-n/2:n/2-1
   acc = acc + im(toUnaliasedIdx(v,n),toUnaliasedIdx(u,n))*dirichlet(x-u,m);
end
z=acc;
