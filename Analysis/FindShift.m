function [ shift ] = FindShift( y,y_ref )
%FINDSHIFT Finds the optimal shift in terms of PSNR between two images
[m,n] = size(y);
[M,N] = size(y_ref);

[~,max_ind] = max(vec(xcorr2(y,y_ref)));
[i,j] = ind2sub([m+M-1,n+N-1],max_ind);


shift = zeros(2,1);
shift(1) = M-i;
shift(2) = N-j;

end

