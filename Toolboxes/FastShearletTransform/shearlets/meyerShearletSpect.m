function Psi = meyerShearletSpect( x, y, a, s, realCoefficients, meyeraux_handle, scaling)
%MEYERSHEARLETSPECT compute the spectrum of the meyer shearlet
% Computes the spectrum of the shearlet "meyerShearlet" for given scale a,
% shear s on the grid spanned by x and y.
% With meyeraux_handle different auxiliary functions can be selected.
%
% INPUT:
%  x	            (meshgrid) the meshgrid for the x-axis
%  y	            (meshgrid) the meshgrid for the y-axis
%  a                (real) scale
%  s                (real) shear
%  realCoefficients (bool) real/complex shearlets
%  meyeraux_handle  (function handle) auxiliary function
%  scaling          ('scaling') compute the respective scaling function
%
% OUTPUT:
%  Psi	(matrix) spectrum
%
%--------------------------------------------------------------------------
% Sören Häuser ~ FFST ~ 2014-07-22 ~ last edited: 2014-07-22 (Sören Häuser)

    %compute scaling function
	if(nargin>6)
		if(strcmp(scaling,'scaling'))
			%cones
			C_hor = abs(x) >= abs(y); %with diag
			C_ver = abs(x) < abs(y);
			Psi = meyerScaling(x,meyeraux_handle) .* C_hor + meyerScaling(y,meyeraux_handle) .* C_ver;
			return;
		end
	end
	
	%set default values
	if(nargin<6)
		meyeraux_handle = @meyeraux;
	end
	if(nargin<5)
		realCoefficients = 1;
	end

	%compute scaling and shearing
	y = s .* sqrt(a) .* x + sqrt(a) .* y;
	x = a .* x;

	%set values with x=0 to 1 (for division)
	xx = (abs(x)==0) + (abs(x)>0).*x;
	
	%compute spectrum
	Psi = meyerWavelet(x,realCoefficients,meyeraux_handle).*bump(y./xx,meyeraux_handle);
end
