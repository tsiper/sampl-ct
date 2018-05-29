function [Y,flag,residual,iter] = fippft3(pp,ErrTol,MaxIts,verbose,precision)
%
% Fast inverse pseudo-polar Fourier transform.
% See ippft3 for more information
%
% The function tries to load a precomputed convolution filter that
% corresponds to the dimensions of the given pseudo-polar array pp.
% If this filter does not exist, the function generates and saves it to the
% disk. The filename of the saved filter is "ppfiltNN.dat". 
% Generating the filter might take a long time, but subsequent calls with
% the same dimensions would be fast.
%
% Yoel Shkolnisky, December 2010.
%
% Revisions:
% Yoel Shkolnisky 19/05/2013 OptimizedprecondAdjPPFT3 was renamed to
%      precondAdjPPFT3. Replaced call accordingly.

if nargin<5
   precision='double';
end

if nargin<4
   verbose = 0;
end

if nargin<3
   MaxIts = 10;
end

if nargin<2
   ErrTol = 1.e-2;
end

n=verifyPP(pp);
filtname=getppfiltname(n,precision);
H=loadppftfilter(filtname);
if ~isstruct(H)
    fprintf('Filter not found. Generating filter\n');
    H=makeppftfilter(n,1,precision);
    saveppftfilter(H,filtname);
end
temp = precondadjppft3(pp);
[Y,flag,residual,iter] = CG('FPtP',temp,{H},ErrTol,MaxIts,zeros(size(temp)),verbose);

if flag
   warning ('Inversion did not converge. Residual error %-2.5e',residual);
end
