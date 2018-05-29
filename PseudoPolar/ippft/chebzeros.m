function t=chebzeros(n)
%
% Compute all zero of the Chebyshev polynomial of order n
%
% Parameters:
%   n   Order of the polynomial. n is an integer. n>=0.
%
% Returned value:
%   t   Array with the n zeros.
%
% Yoel Shkolnisky 22/09/04

if (n<0)
    error('n should be a positive integer');
end

t=-cos((2*[1:n]-1)/n*(pi/2));
