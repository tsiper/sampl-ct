function u=topmul(D,v);
%
% Mulitply a toeplitz matrix, whose diagonal form is D, with the vector v,
% in  O(n log n) operations. 
%
% Input parameters
%   D   Diagonal form of the toeplitz matrix to multiply.
%   v   The vector to multiply.
%
% Output parameters
%   u   The product of the toeplitz matrix that corresponds to D with the
%       vector v. 
%
% Yoel Shkolnisky 04/10/04

n=length(v);
v=v(:);

if length(D)~=2*n
    error('D must be 2*length(v)')
end

v = [v ; zeros(n,1)];
fv = fft(v);
u=ifft(D.*fv);
u=u(1:n);