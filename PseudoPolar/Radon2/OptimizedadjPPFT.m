% function im = OptimizedAdjPPFT(pp1,pp2)
%
% Optimized version of adjppft, the adjoint operator of the
% pseudo-polar Fourier transform operator.
%
% See adjppft.m for more information.
%
% Yoel Shkolnisky 22/10/01

function im = OptimizedAdjPPFT(pp1,pp2)

%Check if the input is of valid size 2x(2n+1)x(n+1)
s1 = size(pp1);
s2 = size(pp2);

if (sum(s1-s2)~=0)
   error('pp1 and pp2 must have the size');
end

if (mod(s1(1),2)~=1) | (mod(s1(2),2)~=1)
   error('pp1 and pp2 must be of size (2n+1)x(n+1)');
end

n = [(s1(1)-1)/2;s1(2)-1];
if (n(1)~=n(2))
   error('Input parameter must be of size (2n+1)x(n+1)');
end

n=n(1);
m=2*n+1;
alpha = 2*(n+1)/(n*m);

%Compute the adjoint of PP1
tmp = zeros(2*n+1,n);
for k=-n:n
   u = fliplr(pp1(toUnaliasedIdx(k,2*n+1),:));
   v = cfrft(u,-k*alpha);
   tmp(toUnaliasedIdx(k,2*n+1),:) = v(1:n);
end

tmp = m*icfftd(tmp,1);
adjpp1 = flipud(tmp(n/2+1:3*n/2,:));

%Compute the adjoint of PP2

tmp = zeros(2*n+1,n);
for k=-n:n
   u = fliplr(pp2(toUnaliasedIdx(k,2*n+1),:));
   v = cfrft(u,-k*alpha);
   tmp(toUnaliasedIdx(k,2*n+1),:) = v(1:n);
end
% To follow the code in adjPPFT we should have transposed each row before we assign it to tmp 
% (and creating an array of size nx(2n+1)). Then we had to apply cfft along the rows.
% To save operations, we assign the vector v to tmp without transpose, apply cfft along columns
% and then transpose the entire matrix at once.

tmp = m*icfftd(tmp,1);
tmp = tmp.';
adjpp2 = flipud(tmp(:,n/2+1:3*n/2,:));

%Combine both adjoints
im = adjpp1+adjpp2;
