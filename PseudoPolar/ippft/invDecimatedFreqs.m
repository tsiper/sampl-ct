function rim=invDecimatedFreqs(fim)
%
% Recover the image whose decimated Fourier transform is fim.
% The input is of size (n+1)x(n+1). 
% The function computed the toeplitz matrix that corresponds to the
% transform that map rim to fim, computes its inverse (O(n^2) operations
% for rim of size nxn), and applies this matrix to all rows and columns of
% fim. Each application requires O(nlogn) operations and hence the total
% complexity for inverting fim is O(n^2logn).
%
% Input parameters
%   fim    Frequnecy image to recover.
%
% Output parameters
%   rim    The recovered image.
% 
% Yoel Shkolnisky 04/10/04


n = size(fim);

if (n(1)~=n(2))
    error('Input matrix must be sqaure');
end

n=n(1)-1;
c=zeros(n,1);
r=zeros(1,n);
for k=-n/2:n/2-1
    for l=-n/2:n/2
        c(k+n/2+1)=c(k+n/2+1)+exp(4*pi*i*l*(-n/2-k)/(2*n+1));
        r(k+n/2+1)=r(k+n/2+1)+exp(4*pi*i*l*(k+n/2)/(2*n+1));
    end
end

[m1,m2,m3,m4]=topinv(c,r);
[D1,D2,D3,D4]=topprepinv(m1,m2,m3,m4);

tmpres=zeros(n+1,n);
rim=zeros(n);

% apply inverse to rows
for k=1:n+1
    v=fim(k,:);
    v=adjF(v);
    u=topinvmul(D1,D2,D3,D4,v);
    tmpres(k,:)=reshape(u(:),1,n);
end

% apply inverse to columns
for k=1:n
    v=tmpres(:,k);
    v=adjF(v);
    u=topinvmul(D1,D2,D3,D4,v);
    rim(:,k)=u;
end