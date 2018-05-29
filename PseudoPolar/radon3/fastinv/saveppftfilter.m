function saveppftfilter(H,fname)
%
% The sturcture of the file is
%     n          INTEGER    Size of the filter is n^3.
%     timestamp  REAL*8     Date and time when the filter was created.
%     precond    INTEGER    nonzero if the kernel includes a
%                           preconditioner.
%     precision  INTEGER    0 for single, 1 for double.
%     data       n^3*REAL*8  Data of the filter. The data is complex so
%                the length is 2xn^3..
%
% Yoel Shkolnisky, December 2010.
%
% Revisions:
% April 2014 Y.S. Add precision variable.

timestamp=now;
n=size(H.filter,1)/2;

if isa(H.filter,'single')
    precision=0;
    precision_fmt='real*4';
elseif isa(H.filter,'double')
    precision=1;
    precision_fmt='real*8';
else
    error('Unknown filter precision');
end

if (2*n~=size(H.filter,2)) || (2*n~=size(H.filter,3))
    error('Filter must be a volume of size 2nx2nx2n');
end

fid=fopen(fname,'w');
fwrite(fid,n,'int');
fwrite(fid,timestamp,'real*8');
fwrite(fid,precision,'int');
fwrite(fid,H.precond,'int');
fwrite(fid,H.filter(:),precision_fmt); 
fclose(fid);