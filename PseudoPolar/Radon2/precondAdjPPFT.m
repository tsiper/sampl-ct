% function im = PrecondAdjPPFT(pp1,pp2) 
%
% Computes the preconditioned adjoint of the pseudo-polar Fourier transform.
%
% pp1,pp2  The pseudo-polar sections.
%          pp1 and pp2 must be of size (2n+1)x(n+1) as results from PPFT.
%
% Returns the image (im) that is the preconditioned adjoint of the PPFT.
% See differences between this function and the function optimizedadjppft.m
% to see how the preconditioner is defined.
% See precond.m for an explicit form of the preconditioner.
%
% Yoel Shkolnisky 17/12/02

function im = precondAdjPPFT(pp1,pp2)

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
   if k==0
      mult = 1/(m^2);
   else
      mult = abs(k*alpha);
   end
   
   u = fliplr(pp1(toUnaliasedIdx(k,2*n+1),:));
   v = mult.*cfrft(u,-k*alpha);
   tmp(toUnaliasedIdx(k,2*n+1),:) = v(1:n);
end

tmp = m*icfftd(tmp,1); %Inverse FFT along columns
adjpp1 = flipud(tmp(n/2+1:3*n/2,:));

%Compute the adjoint of PP2

tmp = zeros(2*n+1,n);
for k=-n:n
   if k==0
      mult = 1/(m^2);
   else
      mult = abs(k*alpha);
   end

   u = fliplr(pp2(toUnaliasedIdx(k,2*n+1),:));
   v = mult.*cfrft(u,-k*alpha);
   tmp(toUnaliasedIdx(k,2*n+1),:) = v(1:n);
end
% To follow the code in adjPPFT we should have transposed each row before we assign it to tmp 
% (and creating an array of size nx(2n+1)). Then we had to apply cfft along the rows.
% To save operations, we assign the vector v to tmp without transpose, apply cfft along columns
% and then transpose the entire matrix at once.

tmp = m*icfftd(tmp,1); %Inverse FFT along columns
tmp = tmp.';
adjpp2 = flipud(tmp(:,n/2+1:3*n/2,:));

%Combine both adjoints
im = (adjpp1+adjpp2)/(m^2);

%Revision record
% 15/1/03   Yoel Shkolnisky     Use cfftd instead of column-wise cfft
