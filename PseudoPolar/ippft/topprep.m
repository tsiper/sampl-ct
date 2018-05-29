function D=topprep(c,r);
%
% Prepare the Toeplitz matrix, whose first column is c and first row is r, 
% for fast multiplication by a vector. The function diagonalizes the
% toeplitz matrix given by c and r
%
% Input parameters
%    c  First column of the toeplitz matrix.
%    r  First row of the toeplitz matrix.
%
% Output parameters
%    D  Diagonal form of the toeplitz matrix. This data used for fast 
%       muliplication of the matrix Toeplitz(c,r) by an arbitrary vector. 
%
% Yoel Shkolnisky 04/10/04


% Algorithm description:
%
% The function embeds the toeplitz matrix, given by c and r, in a circulat
% matrix, diagonalizes it, and returns the eigenvalus of the circulat
% matrix.

if (floateq(c(1),r(1)))
    error('First element of column does not match first element of row');
end

c=c(:);
r=r(:);
n=length(c);
m=length(r);

circmat=[c; 0; r(n:-1:2)];
D=fft(circmat);


