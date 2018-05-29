function testippft
%
% Test the function ippt.
% The function generates a matrix, computes its pseud-polar Fourier
% tarnsform, inverts it, and compares the results.
%
% The output of the function should be very close to zero.
%
% Yoel Shkolnisky 11/10/04

EPS=1.0e-10;
n=128;
im=magic(n);

[pp1,pp2]=PPFT(im);
rim=ippft(pp1,pp2,EPS);
[v,I]=max(abs(im(:)-rim(:)));
v/abs(im(I))