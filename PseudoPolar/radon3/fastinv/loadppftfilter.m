function H=loadppftfilter(fname,verbose)
%
% Load precomputed filter used for fast inversion of the pseudo-polar
% Fourier transform.
%
% Set verbose to non-zero to print filter header.
% On error (for example, if the file does not exist), the functions returns
% -1.
%
% See saveppftfilter for more informations.
%
% Yoel Shkolnisky, December 2010.

if ~exist('verbose','var')
    verbose=0;
end

fid=fopen(fname,'r');

if fid<0
    H=-1;
    return;
end;

n=fread(fid,1,'int');
timestamp=fread(fid,1,'real*8');
precision_flag=fread(fid,1,'int');
precond=fread(fid,1,'int');

if precision_flag==0
    precision='single';
    precision_fmt='real*4';
elseif precision_flag==1
    precision='double';
    precision_fmt='real*8';
else
    error('Unknown precision');
end

if verbose
    fprintf('n=%d\n',n);
    fprintf('Filter created on %s\n',datestr(timestamp));
    fprintf('Precision = %s\n',precision);
    fprintf('precond = %d\n',precond);
end

h=fread(fid,(2*n)^3,precision_fmt);

fclose(fid);

if precision_flag==0
    h=single(h);
end

H.precond=precond;
H.filter=reshape(h,2*n,2*n,2*n);

