function Psi = meyerSmoothShearletSpect( x, y, a, s, realCoefficients, shearletArg, scaling)
%MEYERSMOOTHSHEARLETSPECT compute the spectrum of the smooth meyer shearlet
% Computes the spectrum of the smooth variant of the shearlet "meyerShearlet" for given scale a,
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
% REFERENCES
%  construction based on ideas by 
%  Kanghui Guo, and Demetrio Labate. 
%  "The construction of smooth Parseval frames of shearlets." 
%   Mathematical Modelling of Natural Phenomena 8.01 (2013): 82-105.
%
%--------------------------------------------------------------------------
% Sören Häuser ~ FFST ~ 2014-07-22 ~ last edited: 2014-07-22 (Sören Häuser)

	if(nargin>6)
		if(strcmp(scaling,'scaling'))
			Psi = meyerScaling(x,shearletArg) .* meyerScaling(y,shearletArg);
			return;
		end
	end
	
	if(~realCoefficients)
		error('Complex shearlets not supported for smooth Meyer shearlets!');
	end

	%compute scaling and shearing
	asy = s .* sqrt(a) .* x + sqrt(a) .* y;
	y = a .* y;
	x = a .* x;

	%set values with x=0 to 1 (for division)
	%xx = (abs(x)==0) + (abs(x)>0).*x;
	
	%compute spectrum
	W = sqrt((meyerScaling(2^(-2)*x,shearletArg).*meyerScaling(2^(-2)*y,shearletArg)).^2 - (meyerScaling(x,shearletArg).*meyerScaling(y,shearletArg)).^2);
	Psi = W .* bump(asy./x,shearletArg);
    
    %reset NaN to 0
    Psi(isnan(Psi)) = 0;
end

