function rim = ippftd(pp1,pp2,EPS);
%
% Inverse pseudo-polar Fourier transform.
%
% Direct algorithm that computes the inverse pseudo-polar Fourier
% transform.
%
% Input parameters:
%   pp1,pp2   Pseudo-polar sectors.
%   EPS       Required accuracy of the computation
%
% Output parameters:
%   rim       recovered image
%
% Yoel Shkolnisky 24/01/05

if nargin<3
   EPS = 1.0E-7;
end

n=checkInput(pp1,pp2);
[c,d]=loadConsts(n);

% Prepare output array
fim = zeros(n+1); 

% Scan the pseudo-polar arrays in an "onion-peel" order and recover the
% Fourier samples of each row/column
for j=n/2:-1:1
% recover rows j and -j from pp1
    fu = pp1(toUnaliasedIdx(-2*j,2*n+1),:);
    u  = fim(toUnaliasedIdx(-j,n+1),:);
    dj=d(1,n/2-j+1,:);
    cl=c(1,n/2-j+1,:);
    r1 = recover(fu,u,-2*j,1,reshape(dj,1,n),reshape(cl,1,n),EPS);
%
    fu = pp1(toUnaliasedIdx(2*j,2*n+1),:);
    u  = fim(toUnaliasedIdx(j,n+1),:);
    dj=d(3,n/2-j+1,:);
    cl=c(3,n/2-j+1,:);
    r2 = recover(fu,u,2*j,3,reshape(dj,1,n),reshape(cl,1,n),EPS);
%
% recover columns j and -j from pp2
    fv = (pp2(toUnaliasedIdx(-2*j,2*n+1),:)).';
    v  = fim(:,toUnaliasedIdx(-j,n+1));
    dj=d(2,n/2-j+1,:);
    cl=c(2,n/2-j+1,:);
    c1 = recover(fv.',v.',-2*j,2,reshape(dj,1,n),reshape(cl,1,n),EPS);
%
    fv = (pp2(toUnaliasedIdx(2*j,2*n+1),:)).';
    v  = fim(:,toUnaliasedIdx(j,n+1));
    dj=d(4,n/2-j+1,:);
    cl=c(4,n/2-j+1,:);
    c2 = recover(fv.',v.',2*j,4,reshape(dj,1,n),reshape(cl,1,n),EPS);    
%
% store the recovered rows and columns in the output array
    fim(toUnaliasedIdx(-j,n+1),1:n)  = r1;
    fim(toUnaliasedIdx(j,n+1),2:n+1) = r2;
    fim(2:n+1,toUnaliasedIdx(-j,n+1))  = c1.';
    fim(1:n,toUnaliasedIdx(j,n+1)) = c2.';
end
%
% recover row 0 (no need to recover column 0).
fim(toUnaliasedIdx(0,n+1),toUnaliasedIdx(0,n+1))=pp1(toUnaliasedIdx(0,2*n+1),toUnaliasedIdx(0,n+1));

% recover the orignal image by inverting fim
rim=invDecimatedFreqs(fim);
rim=flipud(rim);


function n=checkInput(pp1,pp2)
% Check that the arrays pp1 and pp2 are properly structured. pp1 and pp2
% should be of size 2x(2n+1)x(n+1). If everything is ok with pp1 and pp2
% the function returns n.
%
s1 = size(pp1);
s2 = size(pp2);

if (sum(s1-s2)~=0)
   error('pp1 and pp2 must have the size');
end

if (mod(s1(1),2)~=1) | (mod(s1(2),2)~=1)
   error('pp1 and pp2 must be of size (2n+1)x(n+1)');
end

n = [(s1(1)-1)/2;s1(2)-1];
if (n(1)~=n(2))
   error('Input parameter must be of size (2n+1)x(n+1)');
end

n=n(1);

%
%
%
%
function [c,d]=loadConsts(n)
%
% Try to load the constants from the disk. If the constants do not exist on
% the disk, compute them using a slow algorithm and print a warning.
%
cname=sprintf('c_%d.mat',n);
dname=sprintf('d_%d.mat',n);

loaded=0;
if exist(cname) & exist(dname)
    load(cname);
    load(dname);
    
    % check that the tables c and d were indeed loaded
    if (exist('c')==1) & (exist('d')==1)
        loaded=1;
    end
end

if ~loaded
    warning('Tables not found on disk. Generating inversion tables...');
    [c,d]=ippftconsts(n);
end

%
%
%
%
function w=recover(fu,u,j,type,dj,cl,EPS);
% Respacing of the frequnecy samples to the cartesian grid.
% The vector fx contains samples fractional samples given by
%        n/2-1
% fu(l) =  sum x(u)exp(-2*pi*i*alpha*u*l/m) 
%        u=-n/2
% where alpha=-2j/n and l=-n/2...n/2.
%
% u contains the samples of the Cartesian Fourier grid that was recovered 
% in previous iterations.
%
% The function returns the vector y given by
%         n/2-1
% w(l) =  sum x(j)exp(-2*pi*i*(2*l)u/m).
%        u=-n/2

n=length(fu)-1;
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
    if (type==1) | (type==3)
        idx2=round([-n/2+N(j/2,n+1)/2+1:n/2-N(j/2,n+1)/2]/(alpha/2));
    else
        idx2=round([-n/2+N(j/2,n+1)/2:n/2-N(j/2,n+1)/2-1]/(alpha/2));
    end
elseif (alpha==2)
    % Handling types 2 and 4 is different than type 1 and 3 because of the
    % flip (reverse angle order) in pp2.
    if (type==1) | (type==3)
        idx2=round([-n/2+N(j/2,n+1)/2:n/2-N(j/2,n+1)/2-1]/(alpha/2));
    else
        idx2=round([-n/2+N(j/2,n+1)/2+1:n/2-N(j/2,n+1)/2]/(alpha/2));
   end
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

% Correspondence test of x and y. 
% Used for validation of the code.
% Remove in the final version.
xs=sort(x);
ys=sort(y);
lx=length(x);
ly=length(y);

if (lx~=ly)
    warning('x and y are not of the same length')
end
if (xs(1)~=ys(1))
    warning('First element in x and y is different')
end
if (xs(lx)~=ys(lx))
    warning('Last element in x and y is different')
end
% END of validation code

f = [u(idx1+n/2+1) fu(idx2+n/2+1) u(idx3+n/2+1)];
w=fastfmmresample(f,y,x,dj,cl,EPS);

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
