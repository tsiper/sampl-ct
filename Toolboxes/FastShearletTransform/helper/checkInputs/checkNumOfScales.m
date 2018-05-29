function TF = checkNumOfScales(x)
%CHECKNUMOFSCALES check NumOfScales
%
% INPUT:
%  x	            numOfScales
%
% OUTPUT:
%  TF	            (bool) result
%
%--------------------------------------------------------------------------
% Sören Häuser ~ FFST ~ 2014-07-22 ~ last edited: 2014-07-22 (Sören Häuser)

   TF = false;
   if(isempty(x))
       TF = true;
   elseif ~(isscalar(x) || ndims(x == 3))
	   error('numOfScales has to be a scalar or a 3-d-matrix!');
   elseif ~isnumeric(x)
       error('numOfScales has to be numeric!');
   else
       TF = true;
   end
end

