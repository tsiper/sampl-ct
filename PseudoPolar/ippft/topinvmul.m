function u=topinvmul(D1,D2,D3,D4,v)
%
% Mulitply the inverse toeplitz matrix, where the diagonal form of its
% Gohberg-Semencul factors are given by D1,D2,D3,and D4, with the
% vector v, in O(n log n) operations. 
% The factors D1,D2,D3,D4 are returned by the function topprepinv.
%
% Input parameters
%  D1 D2 D3 D4     Diagonal form of the Gohberg-Semencul factors of the
%                  inverse topelitz matrix.
% v                Vector to multiply.
%
% Output parameters
% u                The result of inv(T)*v where T is inv(T) is the inverse
%                  toeplitz matrix.
%
% Yoel Shkolnisky 04/10/04

u1=topmul(D2,v);
u1=topmul(D1,u1);

u2=topmul(D4,v);
u2=topmul(D3,u2);

u=u1-u2;
