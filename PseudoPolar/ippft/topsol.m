function x=topsol(c,r,y)
%
% Solve a toeplitz system.
% The toeplitz matrix is given by the column c and the row r.
% The function solves the system 
%     Toeplitz(c,r)*x=y
% in O(n^2) operations using the Levinson iterations.
%
% Input parameters
%    c  First column of the toepliz matrix.
%    r  First row of the toepliz matrix.
%
% Output parameters
%    x  The solution of the system Toeplitz(c,r)*x=y
%
% Yoel Shkolnisky 04/10/04

if (floateq(c(1),r(1)))
    error('First element of column does not match first element of row');
end

if (length(c)~=length(r))
    error('c and r must have the same length')
end

if (length(r)==0)
    y=[];
    return
end

% store c and r in a single array R such that R(n+j) is equal to element
% R(j) where
%    +-                                            -+
%    |R(0)   R(-1)   R(-2) ... R(-(N-2))  R(-(N-1)) |
% R= |R(1)   R(0)    R(-1)     R(-(N-3))  R(-(N-2)) |
%    |...                ....                  .... |
%    |R(N-1) R(N-2) R(N-3) ... R(1)       R(0)      |  
%    +-                                            -+
n=length(c);
c=c(:);
r=r(:);
r=[flipud(r) ; c(2:n)];

% Initialize iterations
x=zeros(n,1);
g=zeros(1,n);
h=zeros(1,n);

x(1)=y(1)/r(n);
if (n==1)
    return
end

g(1)=r(n-1)/r(n);
h(1)=r(n+1)/r(n);

% Iterations
for m=1:n
    m1=m+1;
    sxn=-y(m1);
    sd=-r(n);
    for j=1:m
        sxn=sxn+r(n+m1-j)*x(j);
        sd=sd+r(n+m1-j)*g(m-j+1);
    end
    if (abs(sd)<1.0e-15)
        error('Singular principal minor');
    end
    
    x(m1)=sxn/sd;
    for j=1:m
        x(j)=x(j)-x(m1)*g(m-j+1);
    end
    
    if (m1==n) 
        return
    end
    
    sgn=-r(n-m1);
    shn=-r(n+m1);
    sgd=-r(n);
    for j=1:m
        sgn=sgn+r(n+j-m1)*g(j);
        shn=shn+r(n+m1-j)*h(j);
        sgd=sgd+r(n+j-m1)*h(m-j+1);
    end
    
    if (abs(sgd)<1.0e-15)
        error('Sinular principal minor');
    end
    
    g(m1)=sgn/sgd;
    h(m1)=shn/sd;
    k=m;
    m2=(m+1)/2;
    pp=g(m1);
    qq=h(m1);
    for j=1:m2
        pt1=g(j);
        pt2=g(k);
        qt1=h(j);
        qt2=h(k);
        g(j)=pt1-pp*qt2;
        g(k)=pt2-pp*qt1;
        h(j)=qt1-qq*pt2;
        h(k)=qt2-qq*pt1;
        k=k-1;
    end
end