function [ x_hat, Params ] = EqualizeImage( x,x_true )
%EQUALIZEIMAGE Finds the optimal transformation to transform x to have the best
%PSNR possible when comparing to x_true
% NOTE: We assume that the images are strictly positive

% We find the best affine transformation for moving x to x_true
a = fminsearch(@(a) -psnr(a(1)*x+a(2),x_true),[1,0]);

% Saving our extracted params
Params.Amplitude = a(1);
Params.Offset    = a(2);

% Applying the transformation
x_hat = a(1)*x+a(2);

% We measure the intensity of the negative terms now
Params.NegativeNorm = norm(x_hat(x_hat<0),'fro')/sqrt(numel(x_hat));

Params.PSNR_preNeg = psnr(x_hat,x_true);

% Now we null out the negative terms and the values above 1
x_hat(x_hat<0) = 0;
x_hat(x_hat>max(x_true(:))) = max(x_true(:));

% Measuring before and after PSNR
Params.PSNR_old = psnr(x,x_true);
Params.PSNR = psnr(x_hat,x_true);
Params.SSIM = ssim(x_hat/max(x_true(:)),x_true/max(x_true(:)));


end

