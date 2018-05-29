function v=optimizedtansum(alpha,x,y,eps)
%
% An optimized version of the function tansum. See tansum for more details.
%
% Input parameters:
%    alpha   Weights in the sum above
%    x       Points at which the weights are given.
%    y       Points where to evaluate the sum v.
%    eps     Required accuracy. Default value 1.0e-8.
%
% Returned value:
%    v       An array of the same length as y where v(k) is the value of
%            the sum above at the point y(k).
%
% Yoel Shkolnisky 21/09/04

if (nargin < 3)
    error('Not enough input parameters')
end

if (nargin < 4)
    eps=1.0e-8;
end

if (length(alpha)~=length(x))
    error('Arrays alpha and x must have the same length');
end

% Find the smallest interval the contains all points x and y
lowpoint=min([x y]);
highpoint=max([x y]);

% Sort the arrays alpha and x
[x,sortidx] = sort(x);
alpha=alpha(sortidx);

% Set sizes and legnths
K=length(y);
N=length(x);
depth=floor(log2(N)/512);                        %Number of levels in the algorithm
depth=max(depth,1);
v=zeros(1,K);
rmin=(highpoint-lowpoint)/2^(depth+1);      %The smallest neigborhood processsed in the algorithm
nmax=ceil(log(6/(eps*(sin(rmin))^2))/log(5)); %Largest number of coeffcients used in any level of the algorithm
n=nmax; % all levels use the same value of n since n does not varies very much
t=chebzeros(n);
tmp=zeros(1,n);

% Precompute once the coefficients used by the polynomial u
ucoefs=zeros(n-1,n);
for k=1:n-1
    ucoefs(k,:)=(cos(k*acos(t(:)))).';
end

for l=2:depth
    % Compute the interpolation coefficients f

    r=(highpoint-lowpoint)/2^(l+1);   %Radius of the interval processed in the current level.
    coefs=zeros(2^l,nmax);
    
    xpointidx=1;
    for i=1:2^l 
        % Center of the current interval. The processed interval is
        % [c-r,c+r].
        c=2*r*(i-1)+lowpoint+r;
           
        % Create the far-field expansion coefficients for the current
        % inteval.
        intervalleft=lowpoint+2*r*(i-1);           %left point of the current interval
        intervalright=intervalleft+2*r;  %right point of the current interval
        
        while (xpointidx<=N) & (x(xpointidx)<=intervalright),
            trueidx=sortidx(xpointidx);
            alphak=alpha(trueidx);
            xk=x(trueidx);            
            coefs(i,1:n)=coefs(i,1:n)+alphak.*(t+3*tan(r)*tan(xk-c))./(3*tan(r)-t*tan(xk-c));           
            xpointidx=xpointidx+1;
        end
    end
    
    % For each point y, sum its iteraction with far intervals that were
    % not processed in previous iterations. For each interval i,these
    % intervals are  exactly intervals i-2 and i+2 at the current level
    % (referred to as interaction list at the original paper)
    
    for k=1:K
        % Find the interval of yk
        interval=floor((y(k)-lowpoint)/(2*r))+1;
        if ((interval-2)>=1)            
            c=2*r*(interval-3)+lowpoint+r; % center of interval i-2
            p=3*tan(r)/(tan(y(k)-c));
            % Compute the polynomial u inline (no procedure call) for optimization
            b=cos([1:n-1]*acos(p));
            tmp=ucoefs(1,:).*b(1);
            for m=2:n-1
                tmp=tmp+ucoefs(m,:).*b(m);
            end
            tmp=2.*tmp/n;
            tmp=tmp+1/n;
            
            v(k)=v(k)+sum(coefs(interval-2,1:n).*tmp);
        end
        if ((interval+2)<=2^l)
            c=2*r*(interval+1)+lowpoint+r; % center of interval i+2
            p=3*tan(r)/(tan(y(k)-c));
            % Compute the polynomial u inline (no procedure call) for optimization
            b=cos([1:n-1]*acos(p));
            tmp=ucoefs(1,:).*b(1);
            for m=2:n-1
                tmp=tmp+ucoefs(m,:).*b(m);
            end
            tmp=2.*tmp/n;
            tmp=tmp+1/n;
            
            v(k)=v(k)+sum(coefs(interval+2,1:n).*tmp);
        end
        
        % Check if we should process also intervals i-3 and i+3
        if ((interval-3)>=1) & (floor((interval-1)/2)-floor((interval-4)/2)==1)
            c=2*r*(interval-4)+lowpoint+r; % center of interval i-3
            p=3*tan(r)/(tan(y(k)-c));
            % Compute the polynomial u inline (no procedure call) for optimization
            b=cos([1:n-1]*acos(p));
            tmp=ucoefs(1,:).*b(1);
            for m=2:n-1
                tmp=tmp+ucoefs(m,:).*b(m);
            end
            tmp=2.*tmp/n;
            tmp=tmp+1/n;       
            
            v(k)=v(k)+sum(coefs(interval-3,1:n).*tmp);
        end

        if ((interval+3)<=2^l) & (floor((interval+2)/2)-floor((interval-1)/2)==1)
            c=2*r*(interval+2)+lowpoint+r; % center of interval i+3
            p=3*tan(r)/(tan(y(k)-c));            
            % Compute the polynomial u inline (no procedure call) for optimization
            b=cos([1:n-1]*acos(p));
            tmp=ucoefs(1,:).*b(1);
            for m=2:n-1
                tmp=tmp+ucoefs(m,:).*b(m);
            end
            tmp=2.*tmp/n;
            tmp=tmp+1/n;   

            v(k)=v(k)+sum(coefs(interval+3,1:n).*tmp);
        end      
    end
end

% Compute which points belong to which interval at the finest level.
% The array finest partition contains for for each of the 2^depth intervals at the finest 
% level the indices of all the x points in that interval.
%
% finestpartition is a 2D array that contains for each interval i
% (i=1..2^depth) the indices of all points x that are contained in this
% interval.

finestpartition = zeros(2^depth,N);
r=(highpoint-lowpoint)/2^(depth+1);
xpointidx=1;

for i=1:2^depth
    % Center of the current interval. The processed interval is
    % [c-r,c+r].
    c=2*r*(i-1)+lowpoint+r;
    listidx=1; %current position in the array the corresponds to the current interval.
           
    intervalleft=lowpoint+2*r*(i-1);           %left point of the current interval
    intervalright=intervalleft+2*r;  %right point of the current interval
        
    % The 1.0e-15 compansates for tiny numerical errors
    % XXX: remove the +1.0e-15, run test13 for n=4 (it won't work) and find
    % another fix for the bug. 
%  while (xpointidx<=N) & (x(xpointidx)<=intervalright+1.0e-15),

   while (xpointidx<=N) & (x(xpointidx)<=intervalright),
        finestpartition(i,listidx)=sortidx(xpointidx);
        listidx=listidx+1;
        xpointidx=xpointidx+1;
    end
end

% Process the last point, which may be unprocessed in the above loop due to
% tiny numerical errors
if (xpointidx==N)
    finestpartition(i,listidx)=sortidx(xpointidx);
end
    

% At the finsest level, compute for each point the interaction with the
% neasrest intervals directly.

for k=1:K
    % Find the interval of yk
	interval=floor((y(k)-lowpoint)/(2*r))+1;
    if interval>2^depth %should happen only due to tiny numerical errors
        interval=2^depth;
    end
    
    % Compute the interaction with the previous interval
	if ((interval-1)>=1)
        j=1;
        while (j<=N) & (finestpartition(interval-1,j)~=0)
            idx=finestpartition(interval-1,j);
            v(k)=v(k)+alpha(idx)/tan(y(k)-x(idx));
            j=j+1;
        end
	end

    %Compute the interaction with the current interval
    j=1;
    while (j<=N) & (finestpartition(interval,j)~=0)
        idx=finestpartition(interval,j);
        v(k)=v(k)+alpha(idx)/tan(y(k)-x(idx));
        j=j+1;
    end
    
    % Compute the interaction with the next interval
    if ((interval+1)<=2^depth)
        j=1;
        while (j<=N) & (finestpartition(interval+1,j)~=0)
            idx=finestpartition(interval+1,j);
            v(k)=v(k)+alpha(idx)/tan(y(k)-x(idx));
            j=j+1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% THIS FUNCTION IS NOT USED - all calls were embeded inline.
function cousins=checkcousins(interval1,interval2,l)
%
% Check if the parents of interval1 and interval 2 are neighbors in level 
% l-1. Both interval1 and interval2 are from level l. Hence interval1 and
% interval2 must be in the range 1..2^l.

% This calls are safer but cost a lot time
%if (interval1<1) | (interval1>2^l) 
%    error('interval1 must be in the range 1...2^l');
%end

% This calls are safer but cost a lot time
%if (interval2<1) | (interval2>2^l) 
%    error('interval2 must be in the range 1...2^l');
%end

% The function assumes that l>0 since l==0 should never be called by the
% calling function.
%

cousins=0;
parentint1=floor((interval1-1)/2)+1;
parentint2=floor((interval2-1)/2)+1;
if (abs(parentint1-parentint2)==1) 
    cousins=1;
end

% THIS FUNCTION IS NOT USED - all calls were embeded inline.
function y=u(n,x,ucoefs)
%
% Compute the polynomial u_{j,n}(x) defined in "Fast Fourier Transform for
% Non-equispaced data II" (Dutt and Rokhlin), ACHA 2, 85-100,1995. The
% function is defined by Eq. 52. u_{j,n}(x) is the j'th Lagrange interpolation
% polynomial at Chebyshev zeros t_{1}...t_{n}. The function computes
% u_{j,n}(x) for j=1:n. ucoefs are the precomputed coefficients use by all
% calls to this function
%
% Input parameters:
%   n       Number of interpolation points.
%   x       Point to evaluate the polynomials.
%   ucoefs  Precomputed coefficients used to compute u.
%
% Output:
%   y     Value of u_{j,n}(x) for j=1:n.
%

y=zeros(1,n);
b=cos([1:n-1]*acos(x));
for k=1:n-1
    y=y+ucoefs(k,:).*b(k);
end
y=2.*y/n;
y=y+1/n;
