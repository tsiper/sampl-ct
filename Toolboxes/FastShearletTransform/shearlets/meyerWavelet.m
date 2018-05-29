function y = meyerWavelet(x,realCoefficients,meyeraux_handle)
%
% compute Meyer Wavelet
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

	if(nargin<3)
		meyeraux_handle = @meyeraux;
	end
	if(nargin<2)
		realCoefficients = 1;
	end

	y = sqrt(abs(meyerHelper(x,realCoefficients,meyeraux_handle)).^2+abs(meyerHelper(2*x,realCoefficients,meyeraux_handle)).^2);

end

%helper function
%--------------------------------------------------------------------------
function y = meyerHelper(x,realCoefficients,shearletArg)

	if(realCoefficients)
		xa = abs(x);
	else
		xa = -x; %consider left and upper part of the image due to first row and column
	end

	int1 = ((xa >= 1) & (xa < 2));
	int2 = ((xa >= 2) & (xa < 4));

	psihat = int1 .* sin(pi/2*shearletArg(xa-1));
	psihat = psihat + int2 .* cos(pi/2*shearletArg(1/2*xa-1));

	y = psihat;
end
%--------------------------------------------------------------------------