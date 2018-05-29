function TF = checkImage(x)
%CHECKIMAGE check image input
%
% INPUT:
%  x	            (matrix) image
%
% OUTPUT:
%  TF	            (bool) result
%
%--------------------------------------------------------------------------
% Sören Häuser ~ FFST ~ 2014-07-22 ~ last edited: 2014-07-22 (Sören Häuser)

   TF = false;
   if ~(ismatrix(x) && min(size(x))>1)
       error('A has to be an image!');
   elseif ~isnumeric(x)
       error('A has to be numeric!');
   elseif isinteger(x)
	   warning('Image is integer valued!');
   else
       TF = true;
   end
end
