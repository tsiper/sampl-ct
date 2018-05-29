% function y=GKN3(x,k)
%
% Application of the operator G(K,N) to the sequence x (of length n).
% The operator GKN in 3-D is defined as follows:
%       1. Apply the inverse Fourier transform to the sequence x.
%       2. Pad the resulting sequence to length 3n+1.
%       3. Apply the fractional Fourier transform with alpha=2k/n.
%       4. Return the n+1 central elements.
%
% x  The sequence to resample using the operator GKN. Can be of odd or even length.
% k  The row whose transform is being computed.
%
% Yoel Shkolnisky 30/01/03

function y=GKN3(x,k)

n=length(x);
if (mod(n,2)==1)
   error('Input sequence must of even length');
end
   
w = icfft(x);
w = [zeros(1,n) w zeros(1,n+1)];
w = cfrft(w,2*k/n);  %optimization - use alpha=2k/(n*m) instead of padding
y = w(n+1:2*n+1); % return n+1 central elements
