%function [res1,res2] = OptimizedPPFT(im);
%
% Optimized algorithm for computing the pseudo-polar Fourier transform.
% The computation requires O(n^2logn) operations and uses further 
% optimizations to reduce the operation count (compared to PPFT.m).
%
% im    The image whose pseudo-polar Fourier transform should be computed.
%       Must of a dyadic square size (2^k x 2^k).
%
% Returns res1 and res2 (of size 2n+1xn+1) that contain the pseudo-polar Fourier 
% transform of the input image im.
% res1 contains the values of PPI1(k,l). res2 contains the values of PPI2(k,l). 
%
% See PPFT.m for more documentation.
%
% Yoel Shkolnisky 22/10/01

function [res1,res2] = OptimizedPPFT(im)

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

% initialize constants and data structures
n=s(1);
m=2*n+1;
alpha = 2*(n+1)/(n*m);
res1 = zeros(m,n+1);
res2 = zeros(m,n+1);

% Part I: Computation of Res1
% padding the y-direction and applying column-wise fft
EI  = [zeros(n/2,n);im;zeros(n/2+1,n)];
FEI = cfftd(EI,1); %FFT along columns

% fractional along rows with alpha=2*k*(n+1)/(n*m). This is equivalent
% to padding u to length 2n+1, applying FRFT with alpha=2*k/n and extracting 
% the n+1 central elements. However, this requires less memory and operations.
% The padding with the single zero is used to generate output of length
% n+1.
for k=-n:n
    u = FEI(toUnaliasedIdx(k,m),:);
      w = cfrft([u 0],k*alpha);
      res1(toUnaliasedIdx(k,m),:) = fliplr(w);
end   

% Part II: Computation of Res2
% padding the x-direction and applying row-wise fft
EI  = [zeros(n,n/2) im zeros(n,n/2+1)];
%FEI = (cfftd((EI).',1)).';
FEI = cfftd(EI,2);

for k=-n:n
    v = FEI(:,toUnaliasedIdx(k,m));
      w = cfrft([v.' 0],k*alpha);
      res2(toUnaliasedIdx(k,m),:) = fliplr(w);
end

%Revision record
% 15/1/03   Yoel Shkolnisky     Used cfftd instead of column-wise cfft
