function [c,d]=ippftconsts(n)
% Prepare the tables with the constants used by the interpolation.
% If the tables exists as file then the functions simply reads them.
% Otherwise, it generates the using a slow algorithm. In such a case a
% warning is issued.
%
% Input parameters
%   n    The size of the the pseudo-polar grid is 2x(2n+1)x(n+1)
%
% Output parameters
%   c    Constants cl. 
%   d    Constants dj. 
% Both arrays c and d are 3D arrays of size 4xn2xn . The vectors d(type,:,:) and
% c(type,:,:), type=1,2,3,4, contain constants which correspond to a
% row/column of the given type, where Type 1 refers to rows from -n/2 to -1; type
% 2 refers to columns from -n/2 to -1; type 3 refers to rows from 1 to n/2;
% and type 4 refers to columns from 1 to n/2. 
%
% Yoel Shkolnisky 11/10/04

d=prepConstD(n);
c=prepConstC(n);

%
%
%
%
function d=prepConstD(n)
d=zeros(4,n/2,n);
for j=n/2:-1:1
    d(1,n/2-j+1,:) = compD(n,-2*j,1);
    d(3,n/2-j+1,:) = compD(n,2*j,3);
    d(2,n/2-j+1,:) = compD(n,-2*j,2);
    d(4,n/2-j+1,:) = compD(n,2*j,4);       
end

%
%
%
%
function dj=compD(n,j,type);
%
% Compute the constant srequired for the interpolation of row/column j in
% the the pseudo-polar grid.
%
m=2*n+1;
alpha=-2*j/n;

%Compute the frequnecies in the vector of length N(k,n+1)+n+1. These
%frequnecies are denoted as the vector y
if (type==1) | (type==4)
    idx1=[-n/2:-n/2+N(j/2,n+1)/2-1];
    idx3=[n/2-N(j/2,n+1)/2+1:n/2-1];
else
    idx1=[-n/2+1:-n/2+N(j/2,n+1)/2-1];
    idx3=[n/2-N(j/2,n+1)/2+1:n/2];
end

if (alpha==-2)
    idx2=round([-n/2+N(j/2,n+1)/2+1:n/2-N(j/2,n+1)/2]/(alpha/2));
elseif (alpha==2)
        idx2=round([-n/2+N(j/2,n+1)/2:n/2-N(j/2,n+1)/2-1]/(alpha/2));
elseif (alpha~=0)
    idx2=round([-n/2+N(j/2,n+1)/2:n/2-N(j/2,n+1)/2]/(alpha/2));
else
    idx2=0;
end

y=[-2*pi*(2)*idx1/m -2*pi*alpha*idx2/m -2*pi*(2)*idx3/m];

y=sort(y);

if (mod(n,2)==1)
    error('The length n must be even');
end

LARGE=1.0E15;
EPS=1.0E-15;
dj=zeros(1,n);

% Compute the factors dj
for j=1:n
    dj(j)=0;
    for k=1:n
        if j~=k
            if abs(y(j)-y(k)) > EPS
                dj(j)=dj(j)-log(sin((y(j)-y(k))/2));
            else
                error('two values of y are too close to each other');
            end
        end
    end
end

dj=exp(dj);

%
%
%
%
function c=prepConstC(n)
c=zeros(4,n/2,n);

for j=n/2:-1:1
    c(1,n/2-j+1,:) = compC(n,-2*j,1);
    c(3,n/2-j+1,:) = compC(n,2*j,3);
    c(2,n/2-j+1,:) = compC(n,-2*j,2);
    c(4,n/2-j+1,:) = compC(n,2*j,4);   
end

%
%
%
%
function cl=compC(n,j,type);
%
% Compute the constant required for the interpolation of row/column j in
% the the pseudo-polar grid.
%

m=2*n+1;
alpha=-2*j/n;

%Compute the frequnecies in the vector of length N(k,n+1)+n+1. These
%frequnecies are denoted as the vector y
if (type==1) | (type==4)
    idx1=[-n/2:-n/2+N(j/2,n+1)/2-1];
    idx3=[n/2-N(j/2,n+1)/2+1:n/2-1];
else
    idx1=[-n/2+1:-n/2+N(j/2,n+1)/2-1];
    idx3=[n/2-N(j/2,n+1)/2+1:n/2];
end

if (alpha==-2)
    idx2=round([-n/2+N(j/2,n+1)/2+1:n/2-N(j/2,n+1)/2]/(alpha/2));
elseif (alpha==2)
        idx2=round([-n/2+N(j/2,n+1)/2:n/2-N(j/2,n+1)/2-1]/(alpha/2));
elseif (alpha~=0)
    idx2=round([-n/2+N(j/2,n+1)/2:n/2-N(j/2,n+1)/2]/(alpha/2));
else
    idx2=0;
end

y=[-2*pi*(2)*idx1/m -2*pi*alpha*idx2/m -2*pi*(2)*idx3/m];

%Compute the frequencies on the Cartesian grid.
if (type==1) | (type==4)
    x=-2*pi*(2)*[-n/2:n/2-1]/m;
else
    x=-2*pi*(2)*[-n/2+1:n/2]/m;
end

if (n~=length(x)) || (n~=length(y))
    error('The arrays f,x,y should be of the same length');
end

if (mod(n,2)==1)
    error('The length n must be even');
end

% sort x and then sort everything to the same order as x.
[x,srtXidx]=sort(x);
[y,srtYidx]=sort(y);

% remove x points that are equal to some y point
idxX=1;
idxY=1;
idxDroppedElems=zeros(1,n); % the indices of the x that appear in y

j=1;
while (idxX<=n) & (idxY<=n)
    if (abs(x(idxX)-y(idxY))<1.0d-15)
        idxDroppedElems(j)=idxX;
        j=j+1;
        idxX=idxX+1;
    elseif x(idxX)<y(idxY)
        idxX=idxX+1;
    else
        idxY=idxY+1;
    end
end

% idxUsed is the set of indices from x which we need to
% compute. In other words, these are the elements that were not dropped in
% the previous iteration.
idxUsed=setdiff(1:n,idxDroppedElems(:));
x=x(idxUsed);
m=length(x);
lendropped=j-1; % number of dropped elements from x

% Now we can resample the polynomial using the modified vector x

LARGE=1.0E15;
EPS=1.0E-15;
clt=zeros(1,m);

% Compute the factors cl
for l=1:m
    %compute cl
    clt(l)=0;
    for k=1:n
        clt(l)=clt(l)+log(sin((x(l)-y(k))/2));
    end
%    
    if (abs(clt(l))>LARGE)
        warning('Value of cl too large...');
    end
end

clt=exp(clt);

% Put zeros in the elements that were not computed by this function.
% These zeros will be dropped when the invserion function gets the table of
% cl. We use zero padding so that all arrays will have the same size,
% although this wastes storage space. 
cl=zeros(1,n);
cl(idxUsed)=clt;

%
%
%
%
function l=N(k,s)
% s is the length of the vector x. This means that s=n+1 where n is even
% and the vector x is indexed by
%   [x(-n/2),...,x(0),...,x(n/2)].
% The parameter k should in the range -n/2<=k<=n/2, or when using the
% argument s: -(s-1)/2 <= k <= (s-1)/2.
l=2*((s-1)/2-abs(k));
