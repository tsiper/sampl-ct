function A = mySquare(A,mx,my,r2,a)
% insert square in image A with central position (mx,my) with side length 2r2. If
% given, values are set to a.
%
% INPUT:
%  A                (matrix) image
%  mx               (int) vertical coordinate (column) for center point
%  my               (int) horizontal coordinate (row) for center point
%  r2               (int) half side length
%  a                (real) value of pixels inside square 
%                          (optional, default 1)
%
% OUTPUT:
%  A				(matrix) image with square
%
%--------------------------------------------------------------------------
% Sören Häuser ~ FFST ~ 2014-07-22 ~ last edited: 2014-07-22 (Sören Häuser)

	if(nargin<5)
		a=1;
	end
	r = fix(r2/2);
	
	ax1 = max(mx-r,1);
	ax2 = min(mx+r,size(A,2));
	
	ay1 = max(my-r,1);
	ay2 = min(my+r,size(A,1));
	
	A(ay1:ay2,ax1:ax2) = a;
end