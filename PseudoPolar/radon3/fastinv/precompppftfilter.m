function precompppftfilter(n,precision)
%
% Compute and save the convolution filter that corresponds to the Gram
% operator of the pseudo-polar Fourier tranform.
%
% The filter is saved in a file named ppfiltnnn_sss.dat, where nnn is the
% size n given as an argument, and sss is 'single' or 'double'.
% Default is double.
%
% Yoel Shkolnisky, December 2010.
%
% Revisions:
% April 2014  Y.S.  Add precision variable.

if nargin<2
    precision='double';
end

filtname=getppfiltname(n,precision);
H=makeppftfilter(n,1,precision);
saveppftfilter(H,filtname);
