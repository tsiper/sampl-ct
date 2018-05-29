function [D1,D2,D3,D4]=topprepinv(m1,m2,m3,m4)
%
% Like topprep but processes 4 toeplitz matrices at once.
%
% Input parameters
%   m1 m2 m3 m4     Cell arrays of size 2 where m_i{1} is the first column
%                   of the toeplitz matrix m_i and m_i{2} is the first row of
%                   the toeplitz matrix m_i.
%
% Output parameters
%  D1 D2 D3 D4      D_i is the diagonal form of the matrix m_i.
%
% Yoel Shkolnisky 04/10/04

m1c=m1{1}; m1r=m1{2};
m2c=m2{1}; m2r=m2{2};
m3c=m3{1}; m3r=m3{2};
m4c=m4{1}; m4r=m4{2};

D1=topprep(m1c,m1r);
D2=topprep(m2c,m2r);
D3=topprep(m3c,m3r);
D4=topprep(m4c,m4r);