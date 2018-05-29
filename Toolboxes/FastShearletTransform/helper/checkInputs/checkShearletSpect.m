function TF = checkShearletSpect(x)
%CHECKSHEARLETSPECT check shearletSpect input
%
% INPUT:
%  x	            (string/fhandle) shearletSpect
%
% OUTPUT:
%  TF	            (bool) result
%
%--------------------------------------------------------------------------
% Sören Häuser ~ FFST ~ 2014-07-22 ~ last edited: 2014-07-22 (Sören Häuser)

	TF = false;
	if ~(ischar(x) || isa(f, 'function_handle'))
		error('shearletSpect has to be a string or a function handle!')
	else
		TF = true;
	end
end

