function [ g ] = adjFam( f,a )
%FA Calculates the adjoint fractional psuedo polar fourier for 
% the function f with paramater a and length m=2*n+1
% The length of f is n+1 and the length of g is n

% % f = f(:);
% % n = length(f)-1;
% % % m = 2*n+1;
% % m = 2*n+1;
% % u = (-n/2:n/2)';
% % g = zeros(m,1);
% % 
% % for k = -n:n
% %     g(k+n+1) = sum( f.*exp(+2j*pi*a*k*u/m) );
% % end
% % 
% % % Cropping back to the n central elements
% % g = g(n/2+1:end-n/2-1);
% % 
% % end
f = f(:);
n = length(f)-1;
m = 2*n+1;

W = exp(+2j*pi*a/m);
A = 1;

% g = chirpZ(f,A,W,n);
g = OptimizedChirpZ(f,A,W,n);

g = g(:);