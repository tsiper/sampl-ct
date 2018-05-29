function [ ] = ShowImage( I ,Title , dbFlag)
%UNTITLED Performs full analysis of the image I
Contrast = 1e-4;
myfft = @(x) db((abs(fftshift(fft2(x)))/numel(x))+Contrast);
if nargin < 2
    Title = '';
end
if nargin == 3
    if strcmpi(dbFlag,'db')
        Eps = 1e-2;
        myfft = @(x) db(abs(fftshift(fft2(x)))/numel(x)+Eps);
    end
end


figure('Name',Title,'Position',[50 50 1280 640]);
gap = [0.1];
marg_h= [0.05 0];
marg_w = [0 0.1];
ax(1) = subtightplot(1,2,1,gap,marg_h,marg_w);imagesc(I);colormap 'bone';freezeColors;
ax(2) = subtightplot(1,2,2,gap,marg_h,marg_w);imagesc(myfft(I));colormap 'default';
if nargin > 1
    suptitle(Title);
end
MyColorbar(ax);

end

