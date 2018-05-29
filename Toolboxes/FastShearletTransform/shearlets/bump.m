function y = bump(x,meyeraux_handle)
% compute the function psi_2^ at given points x
%
% INPUT:
%  x                (vector) grid points
%  meyeraux_handle  (function handle) auxiliary function
%
% OUTPUT:
%  y				(vector) values at given points x
%
%--------------------------------------------------------------------------
% Sören Häuser ~ FFST ~ 2014-07-22 ~ last edited: 2014-07-22 (Sören Häuser)

	if(nargin<2)
		meyeraux_handle = @meyeraux;
	end

	y = meyerBump(1+x,meyeraux_handle).*(x<=0) + meyerBump(1-x,meyeraux_handle).*(x>0);
	y = sqrt(y);

end

%--------------------------------------------------------------------------
function y = meyerBump(x,shearletArg)

	int1 = shearletArg(x).*(x>=0).*(x<=1);
	y = int1 + (x>1);

end
%--------------------------------------------------------------------------