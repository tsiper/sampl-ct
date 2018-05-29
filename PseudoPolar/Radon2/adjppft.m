% function im = adjPPFT(pp1,pp2) 
%
% Computes the adjoint of the pseudo-polar Fourier transform.
%
% pp1,pp2   The pseudo-polar sections resulted from the PPFT (pseudo-polar Fourier transform).
%           pp1 and pp2 must be of size (2n+1)x(n+1) as results from PPFT.
%
% See also PPFT.
%
% Yoel Shkolnisky 9/2/02


function im = adjPPFT(pp1,pp2)

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
adjpp1 = zeros(n);
adjpp2 = zeros(n);

%Compute the adjoint of PP1
tmp = zeros(2*n+1,n);
for k=-n:n
   u = fliplr(pp1(toUnaliasedIdx(k,2*n+1),:));
   v = adjGKN(u,k);
   tmp(toUnaliasedIdx(k,2*n+1),:) = v;
end
tmp = prod(size(tmp))*icfft2(tmp);
adjpp1 = flipud(tmp(n/2+1:3*n/2,:));

%Compute the adjoint of PP2
tmp = zeros(n,2*n+1);
for k=-n:n
   u = fliplr(pp2(toUnaliasedIdx(k,2*n+1),:));
   v = adjGKN(u,k);
   tmp(:,toUnaliasedIdx(k,2*n+1)) = v.';
end
tmp = prod(size(tmp))*icfft2(tmp);
adjpp2 = flipud(tmp(:,n/2+1:3*n/2,:));

%Combine both adjoints
im = adjpp1+adjpp2;