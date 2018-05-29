function y = meyerScaling(x,meyeraux_handle)
%MEYERSCALING compute Meyer scaling function
% Compute mother scaling function for meyer shearlet
%
% INPUT:
%  x                (vector) grid points
%  meyeraux_handle  (function handle) auxiliary function
%
% OUTPUT:
%  y				(vector) values at given points x
%
%--------------------------------------------------------------------------
% Sören Häuser ~ FFST ~ 2014-07-22 ~ last edited: 2014-07-22 (Sören Häuser)

	if(nargin<2)
		meyeraux_handle = @meyeraux;
	end

	xa = abs(x);

	% Compute support of Fourier transform of phi.
	int1 = ((xa < 1/2));
	int2 = ((xa >= 1/2) & (xa < 1));

	% Compute Fourier transform of phi.
	phihat = int1 .* ones(size(x));
	phihat = phihat + int2.* cos(pi/2*meyeraux_handle(2*xa-1));

	y = phihat;
end 