function TF = checkLength(x)
%CHECKLENGTH check coefficient input
%
% INPUT:
%  x	            (vector) size
%
% OUTPUT:
%  TF	            (bool) result
%
%--------------------------------------------------------------------------
% Sören Häuser ~ FFST ~ 2014-07-22 ~ last edited: 2014-07-22 (Sören Häuser)

   TF = false;
   if ~(vector(x) && length(x)==2)
       error('l has to be a vector [rows columns]!');
   elseif ~isnumeric(x)
       error('l has to be numeric!');
   else
       TF = true;
   end
end
