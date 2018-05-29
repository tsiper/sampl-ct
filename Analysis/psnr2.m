function [ psnr_val ] = psnr2( x,y )
%PSNR2 Calculates the best PSNR, even when shifts are required

xycorr_max   = max(vec(xcorr2(x,y)));
DynamicRange = max(y(:)) - min(y(:));
MSE          = (norm(x(:))^2+norm(y(:))^2 - 2*xycorr_max)/numel(y);

if MSE > 0
    psnr_val = 10*log10(DynamicRange^2 / MSE );
else
    psnr_val = inf;
end


end

