function [m1,m2,m3,m4]=topinv(c,r)
%
% factorize the inverse toeplitz matrix using the Gohberg-Semencul
% factorization. This requires O(n^2) operations. 
%
% The function returns four toeplitz matrices which are the
% Gohberg-Semencul factorization.
%
% Input parameters
%    c  First column of the toeplitz matrix.
%    r  First row of the toeplitz matrix.
%
% Output parameters
%    m1 m2 m3 m4     Four toeplitz matrics, which correpond to the
%                    Gohberg-Semencul factorization of inv(Toeplitz(c,r)).
%                    There are cell arrays of size 2 where m_i{1} is the first column 
%                    of the toeplitz matrix m_i and m_i{2} is the first row of
%                    the toeplitz matrix m_i.
%
% Yoel Shkolnisky 04/10/04

if (floateq(c(1),r(1)))
    error('First element of column does not match first element of row');
end
n=length(c);
m=length(r);

if (m~=n)
    error('c and r must have the same length');
end

c=c(:);
r=r(:);

e0=zeros(m,1);
en=zeros(m,1);
e0(1)=1;
en(n)=1;

x=topsol(c,r,e0);
y=topsol(c,r,en);

m1={reshape(x,n,1)./x(1), [1 zeros(1,n-1)]};
m2={[y(n); zeros(n-1,1)], reshape(y(n:-1:1),1,n)};
m3={[0; reshape(y(1:n-1),n-1,1)]./x(1), zeros(1,n)};
m4={zeros(n,1) , [0 reshape(x(n:-1:2),1,n-1)]};

