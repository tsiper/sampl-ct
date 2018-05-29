function [ g ] = Fam( c,a )
%FA Calculates the fractional psuedo polar fourier for 
% the function f with paramater a and length n. m=2*n+1
% The length of g is n+1
% % 
% % c = c(:);
% % n = length(c);
% % m = 2*n+1;
% % 
% % % Padding to length m=2*n+1
% % % c = padarray(c,[n 0]);
% % % m = 2*n;
% % 
% % u = (-n/2:n/2-1)';
% % g = zeros(m,1);
% % 
% % for k = -n:n
% %     g(k+n+1) = sum( c.*exp(-2j*pi*a*k*u/m) );
% % end
% % 
% % % Cropping back to the n+1 central elements
% % g = g(n/2+1:end-n/2);
% % 
% % end
c = c(:);
n = length(c);
m = 2*n+1;

W = exp(-2j*pi*a/m);
A = 1;

% g = chirpZ(c,A,W,n+1);
g = OptimizedChirpZ(c,A,W,n+1);

g = g(:);
