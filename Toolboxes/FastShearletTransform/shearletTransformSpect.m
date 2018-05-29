function [ST,Psi] = shearletTransformSpect(A,varargin)
%SHEARLETTRANSFORMSPECT compute shearlet transform
% Compute the shearlet transform of a given image A. The number of scales
% and a boolean indicating real or complex coefficients are optional
% parameters.
% Using a parameter value list further details can be provided.	
% ST contains the shearlet coefficients in a 3-d-matrix
% where the third index indicates the respective shear (1 for lowpass,
% 2:end for the different shears and scales). The images are ordered
% ascending with the scale and within each scale counter-clockwise with the
% direction of the shear. 
% Psi contains the respective shearlets (in Fourier domain).
%
% INPUT:
%  A                (matrix) image (or data) to transform
%  numOfScales		(int) number of scales OR
%	 				(3-d-matrix) precomputed Psi (optional)
%  realCoefficients (bool) real/complex shearlets  (optional)
%
% OUTPUT:
%  ST				(3-d-matrix) shearlet coefficients
%  Psi				(3-d-matrix) spectrum of shearlets
%
% PARAMETERS: (as optional parameter value list, arbitrary order)
%  'shearletSpect'  (string or function handle) shearlet spectrum
%  'shearletArg'	(arbitrary) further parameters for shearlet
%  'realReal'       (bool) guarantees really real shearlets
%  'maxScale'       ('max','min') maximal or minimal finest scale
%
% EXAMPLES:
%  %compute shearlet transform of image A with default parameters	
%  [ST,Psi] = shearletTransformSpect(A);
%  %this is equivalent to
%  [ST,Psi] = shearletTransformSpect(A,floor(0.5*log2(max(size(A)))),1,'shearletSpect',@meyerShearletSpect,'shearletArg',@meyeraux,'realReal',1,'maxScale','max');
%
%  %compute shearlet transform of image A with precomputed shearlet spectrum
%  [ST,Psi] = shearletTransformSpect(A,Psi);
%
%  %compute sharlet transform of image A with specified number of scales
%  [ST,Psi] = shearletTransformSpect(A,4);
%
%  %compute shearlet transform of image A with own shearlet
%  [ST,Psi] = shearletTransformSpect(A,'shearletSpect',@yourShearletSpect);
%
%--------------------------------------------------------------------------
% Sören Häuser ~ FFST ~ 2014-07-22 ~ last edited: 2014-07-22 (Sören Häuser)

    %% parse input
    p = inputParser;
	
	addRequired(p,'A',@checkImage);
 	addOptional(p,'numOfScales',defaultNumberOfScales(size(A)),@checkNumOfScales);
	addOptional(p,'realCoefficients',1,@isnumeric);
	p = parseShearletParameterInputs( p );	
    
	parse(p,A,varargin{:});
    pp = p.Results;
	
	if(ischar(pp.shearletSpect))
		pp.shearletSpect = str2func(pp.shearletSpect);
    end
    
    %add size of image to struct
    pp.l = size(A);

	%% compute spectra
	if( ndims(pp.numOfScales) == 3 )
		Psi = pp.numOfScales;
	else
		Psi = scalesShearsAndSpectra ( size(A), pp );
    end
    
	%% shearlet transform
	uST = Psi .* repmat(fftshift(fft2(A)),[1,1,size(Psi,3)]);
	ST = ifft2(ifftshift(ifftshift(uST,1),2));
	clear uST;
	
	%due to round-off errors the imaginary part is not zero but very small
	%(~ 10^-18) -> neglect it
	if(pp.realCoefficients && pp.realReal)
  		ST = real(ST);
	end
		
end