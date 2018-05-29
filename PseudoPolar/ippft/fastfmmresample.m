function g=fastfmmresample(f,y,x,dj,cl,EPS)
%
% Resample a trigonometric polynomial from the points y to the points x.
%
% Let f(w) be the trigonometric polynomial given by
%
%       n/2-1
% T(w) = sum a(k)*exp(i*k*w).
%       k=-n/2
% The array f is the sample of T(w) at the points y. The function computes
% the values of T(w) at the points x.
%
% Input parameters:
%   f      Values of T(x) at the points y.
%   y      Points where f is given.
%   x      Points where to resample T(w).
%   cl,dj  Precomputed interpolation constants.
%
% Output parameters:
%   g     Values of T(w) at the points x.
%
% The arrays f,x,y should be of the same even length.
%
% Yoel Shkolnisky  03/10/2004

n=length(f);

if (nargin<5)
    EPS=10.e-8
end

if (nargin<5)
    error('Interpolation constants cl and dj must be supplied')
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
f=f(srtYidx);

% remove x points that are equal to some y point
idxX=1;
idxY=1;
idxDroppedElems=zeros(1,n); % the indices of the x that appear in y
idxEqualTo=zeros(1,n); % the index in the y array of the dropped elements

j=1;
while (idxX<=n) & (idxY<=n)
    if (abs(x(idxX)-y(idxY))<1.0d-15)
        idxDroppedElems(j)=idxX;
        idxEqualTo(j)=idxY;
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
cl=cl(idxUsed);

% Now we can resample the polynomial using the modified vector x

LARGE=1.0E15;
g=zeros(1,n);

% Resample the polynomial
gg=zeros(1,m);
if length(x)>0
	fd=f.*dj;
	sfd=sum(fd);
	b=optimizedtansum(fd,y./2,x./2,EPS);
	for k=1:m
        b(k)=b(k)-i*sfd;
	end
	gg=cl.*b;
end

% For each dropped x, the value of the polynomial is simply the
% corresponding value from f. 
% Create the output array by appending the computed values to the relevant
% values from f.
g(idxUsed)=gg;
for k=1:lendropped
    g(idxDroppedElems(k))=f(idxEqualTo(k));
end

% reorder the output array to the original order.
g=g(srtXidx);
