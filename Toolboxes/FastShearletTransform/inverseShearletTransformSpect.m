function A = inverseShearletTransformSpect (ST, varargin)
%INVERSESHEARLETTRANSFORMSPECT compute inverse shearlet transform
% Compute the inverse shearlet transform for given shearlet coefficients ST.
% If the shearlet spectra are not given they are computed using parameters 
% guessed from the coefficients.
% The parameters 'shearletSpect', 'shearletArg' and 'maxScale' cannot be 
% guessed and have to be provided if not the default ones.
%
% INPUT:
%  ST				(3-d-matrix) shearlet transform
%  Psi				(3-d-matrix) spectrum of shearlets (optional)
%
% OUTPUT:
%  A				(matrix) reconstructed image
%
% PARAMETERS: (as optional parameter value list, arbitrary order)
%  'shearletSpect'  (string or function handle) shearlet spectrum
%  'shearletArg'	(arbitrary) further parameters for shearlet
%  'maxScale'       ('max','min') maximal or minimal finest scale
%
%--------------------------------------------------------------------------
% Sören Häuser ~ FFST ~ 2014-07-22 ~ last edited: 2014-07-22 (Sören Häuser)

    %% parse input
    p = inputParser;

	addRequired(p,'ST',@checkCoefficients);
	addOptional(p,'Psi',NaN); %will be computed later
	p = parseShearletParameterInputs( p );

	parse(p,ST,varargin{:});
    pp = p.Results;
	
	if (isnan(pp.Psi))
		%numOfScales
		%possible: 1, 4, 8, 16, 32,...
		% -> -1 for lowpass
		% -> divide by for (1, 2, 4, 8, ...
		% -> +1 results in a 2^# number -> log returns #
		pp.numOfScales = log2((size(ST,3)-1)/4 + 1);
		
		%realCoefficients
		pp.realCoefficients = isreal(ST(:,:,2));
		
		%realReal
		pp.realReal = isreal(ST(:,:,end));
		
		%compute spectra
		pp.Psi = scalesShearsAndSpectra ([size(ST,1),size(ST,2)],pp);
	end

	%% inverse shearlet transform
	A = sum(fftshift(fftshift(fft2(ST),1),2).*pp.Psi,3);
	A = real(ifft2(ifftshift(A)));

end