function TF = checkCoefficients(x)
%CHECKCOEFFICIENTS check coefficient input
%
% INPUT:
%  x	            (matrix) coefficients
%
% OUTPUT:
%  TF	            (bool) result
%
%--------------------------------------------------------------------------
% Sören Häuser ~ FFST ~ 2014-07-22 ~ last edited: 2014-07-22 (Sören Häuser)

   TF = false;
   if ~(ndims(x) == 3)
       error('ST has to be a 3-d matrix of coefficients!');
   elseif ~isnumeric(x)
       error('ST has to be numeric!');
   else
       TF = true;
   end
end
