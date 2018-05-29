% function [re1,res2] = PPFT(im);
%
% Fast algorithm for computing the pseudo-polar Fourier transform.
% The computation requires O(n^2logn) operations.
%
% im    The image whose pseudo-polar Fourier transform should be computed.
%       Must of a dyadic square size (2^k x 2^k).
%
% Returns res1 and res2 (of size 2n+1xn+1) that contain the pseudo-polar Fourier 
% transform of the input image im.
% res1 contains the values of PPI1(k,l). res2 contains the values of PPI2(k,l). 
% The first argument of Res1 and Res2 corresponds to pseudo-radius and the second argument 
% corresponds to pseudo-angle.
% Angle ordering:
%       res1(k,l) - pseudo-polar fourier transform which corresponds to radius "k" and angle
%                   arctan(2l/n). "l" goes from -n/2 to n/2, and therefore, for a constant "k",
%                   res1 corresponds to angles from -pi/4 to pi/4.
%       res2(k,l) - pseudo-polar fourier transform which corresponds to radius "k" and angle
%                   pi/2-arctan(2l/n). "l" goes from -n/2 to n/2, and therefore, for a constant "k",
%                   res1 corresponds to angles from 3pi/4 to pi/4.
% To obtains a continuous change in angle from res1 to res2, one must flip res2 such that it
% corresponds to angles from pi/4 to 3pi/4 (and therefore continuous the angles in res1).
%
%     3pi/4   pi/4
%        \    /
%         \  /  ^
%          \/   |
%          /\   |  angle change in res1 
%         /  \  |
%        /    \
%             -pi/4
%
%
%
%       -------> angle change in res2
%     3pi/4   pi/4
%        \    /
%         \  /  
%          \/   
%          /\    
%         /  \  
%        /    \
%             -pi/4
%
% See thesis' final PDF for more information.
%
% Yoel Shkolnisky 22/10/01

function [res1,res2] = PPFT(im);

% Matlab origin if top-left corner. y-axes goes up-down (values increase as we move down). 
% Mathematical axes goes bottom-up (negative values are at the bottom). We must flip the matrix
% to convert Matlab axes into the mathematical axes so the matrix will match the mathematical equations.
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
n=s(1);  % At this point, all validity test passed and therefore the input image is square. Therefore we can
            % choose n to be the size of either of the dimensions.
m=2*n+1;
res1 = zeros(m,n+1);
res2 = zeros(m,n+1);

% Computation of Res1
EI  = [zeros(n/2,n);im;zeros(n/2+1,n)];
FEI = cfft2(EI);
      
for k=-n:n
    u = FEI(toUnaliasedIdx(k,m),:);
      w = GKN(u,k);
      res1(toUnaliasedIdx(k,m),:) = fliplr(w); %See thesis paper for explanation of the flip
end   


% Computation of Res2
EI  = [zeros(n,n/2) im zeros(n,n/2+1)];
FEI = cfft2(EI);

for k=-n:n
    v = FEI(:,toUnaliasedIdx(k,m));
      w = GKN(v.',k);
      res2(toUnaliasedIdx(k,m),:) = fliplr(w); %See thesis paper for explanation of the flip
end

%Revision record
% 9/12/02   Yoel Shkolnisky     Added comments
