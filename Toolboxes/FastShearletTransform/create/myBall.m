function A = myBall(A,mx,my,r,a)
% insert circle in image A with central position (mx,my) with radius r (2-ball). If
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
		a = 1;
	end

	if(length(r)==1)
		for(i=1:size(A,2))
			for(j=1:size(A,1))
				if(((i-mx)^2+(j-my)^2)<=r^2)
					A(j,i) = a;
				end
			end
		end
	else
		for(i=1:size(A,2))
			for(j=1:size(A,1))
				if((((i-mx)^2+(j-my)^2)<=r(2)^2) && (((i-mx)^2+(j-my)^2)>=r(1)^2))
					A(j,i) = a;
				end
			end
		end
	end
end