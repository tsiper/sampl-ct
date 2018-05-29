function [ psnr_val, snr_val ] = psnr( y1, y2 )
%PSNR Computes the PSNR between two images
% n is the length of one row/column unless 

% Vectorizing the images
y1 = y1(:);
y2 = y2(:);

if size(y1) ~= size(y2), error('psnr - Images aren''t the same size'); end
%     
% if nargin < 3
%     [n,m] = size(y2);
% end
    
DynamicRange = max(y2)-min(y2);
NoiseEnergy  = norm(y1-y2,'fro')^2;
MSE          = NoiseEnergy/numel(y2);
psnr_val     = 10*log10(max(DynamicRange^2,1) / MSE);
snr_val      = 10*log10(norm(y2(:))^2 / NoiseEnergy);

% Because we work on sparse vector we add the full function
psnr_val = full(psnr_val);

end