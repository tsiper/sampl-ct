function p = parseShearletParameterInputs( p )
%PARSESHEARLETINPUTS general input parser
% parse input arguments and set default values
%
% INPUT:
%  p	            (object) inputParser object
%
% OUTPUT:
%  p	            (object) inputParser object with parameter-value pairs
%
%--------------------------------------------------------------------------
% Sören Häuser ~ FFST ~ 2014-07-22 ~ last edited: 2014-07-22 (Sören Häuser)

	validMaxScales = {'min','max'};
	
	addParamValue(p,'shearletSpect',@meyerShearletSpect,@checkShearletSpect);
	addParamValue(p,'shearletArg',@meyeraux);
	addParamValue(p,'realReal',1,@(x) isnumeric(x) || islogical(x));
	addParamValue(p,'maxScale','max',@(x) any(validatestring(x,validMaxScales)));

end

