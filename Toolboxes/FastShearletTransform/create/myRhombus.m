function A = myRhombus(A, mx,my,r, a)
% insert rhombus in image A with central position (mx,my) with radius r (1-ball). If
% given, values are set to a.
%
% INPUT:
%  A                (matrix) image
%  mx               (int) vertical coordinate (column) for center point
%  my               (int) horizontal coordinate (row) for center point
%  r                (int) radius
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

	for(i=1:size(A,2))
		for(j=1:size(A,1))
			if((abs(i-mx)+abs(j-my))<=r)
				A(j,i) = a;
			end
		end
	end
	
end