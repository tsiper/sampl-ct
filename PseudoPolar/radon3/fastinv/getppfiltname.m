function filtname=getppfiltname(n,precision)
%
% Genetrate the name of the file that stores the filter that corresponds to
% size n and given precision.
% 
% Yoel Shkolnisky, December 2010.
%
% Revisions:
% April 2014 Y.S.   Add precision variable.

if nargin<2
    precision='double';
end

filtname=sprintf('ppfilt%d_%s.dat',n,precision);