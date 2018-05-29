% function y=GKN(x,k)
%
% Application of the operator G(K,N) to the sequence x (of length n).
% The operator GKN is defined as follows:
%       1. Apply the inverse Fourier transform to the sequence x.
%       2. Pad the resulting sequence to length 2n+1.
%       3. Apply the fractional Fourier transform with alpha=2k/n.
%       4. Return the n+1 central elements.
%
% x   The sequence to resample using the operator GKN. Can be of odd or even
%     length. Must be a 1-D row vector.
% k   The row whose transform is being computed.
%
% Returns the result of the application of GKN to the sequence x.
% See thesis' final version for more information.
% 
% See Also adjGKN.
%
% Yoel Shkolnisky 22/10/01

function y=GKN(x,k)

n=length(x);
if (mod(n,2)==1)
   error('Input sequence must of even length');
end
   
w = icfft(x);
w = [zeros(1,n/2) w zeros(1,n/2+1)];
w = cfrft(w,2*k/n);  %optimization - use alpha=2k/m instead of padding
y = w(n+1-n/2:n+1+n/2);
