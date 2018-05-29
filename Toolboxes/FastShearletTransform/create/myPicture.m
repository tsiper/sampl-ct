function A = myPicture()
% create image with geometric shapes
%
% INPUT:
%
% OUTPUT:
%  A				(matrix) image with some geometric shapes
%
%--------------------------------------------------------------------------
% Sören Häuser ~ FFST ~ 2014-07-22 ~ last edited: 2014-07-22 (Sören Häuser)

	A = zeros(512,512);
	A = myRhombus(A,348,348,48);
	A = mySquare(A,256-16,348,64);
	A = myBall(A,256+16,348+32,64);
	A = myBall(A,348,128,64);
	A = mySquare(A,128,128,128);
end
