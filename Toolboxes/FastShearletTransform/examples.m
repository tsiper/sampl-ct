% show different examples of computation of shearlet coefficients
% just run the whole script or section by section
%
%--------------------------------------------------------------------------
% Sören Häuser ~ FFST ~ 2014-07-22 ~ last edited: 2014-07-22 (Sören Häuser)

%init
clear all
close all

%load example image
X = double(imread('cameraman.tif'))/255;
% note that conversion to double is necessary, range of image can be
% arbitrary

%% compute shearlet transform
[ST,Psi] = shearletTransformSpect(X);

%% shearlet coefficients and spectra using subplot
idx = 14;
figure
subplot(2,2,1)
h = imagesc(X);
axis image off
colormap(gray)
colorbar
title('original image')

subplot(2,2,2)
imagesc(ST(:,:,idx))
axis image off
colorbar
title('shearlet coefficients')

subplot(2,2,3)
imagesc(Psi(:,:,idx))
axis image off
colorbar
title('shearlet in Fourier domain')

subplot(2,2,4)
imagesc(fftshift(ifft2(ifftshift(Psi(:,:,idx)))))
axis image off
colorbar
title('sherlet in time domain')

%% shearlet coefficients and spectra using pause
figure
for(i=1:size(ST,3))
	imagesc(ST(:,:,i))
	axis image off
	colorbar
    title(sprintf('Shearlet coefficients for idx %d',i))
    disp('press any key to continue...')
	pause
end

figure
for(i=1:size(Psi,3))
	imagesc(Psi(:,:,i))
	axis image off
	colorbar
    title(sprintf('Shearlet in Fourier domain for idx %d',i))
    disp('press any key to continue...')
	pause
end

figure
for(i=1:size(Psi,3))
	imagesc(fftshift(ifft2(ifftshift(Psi(:,:,i)))))
	axis image off
	colorbar
    title(sprintf('Shearlet in time domain for idx %d',i))
    disp('press any key to continue...')
	pause
end

%% show frame tightness and exactness
figure
imagesc(1 - sum(Psi.^2,3))
axis image off
colorbar
title('Frame Tightness')

figure
XX = inverseShearletTransformSpect(ST,Psi);
imagesc(X-XX)
axis image off
colorbar
title('Transform Exactness')

%% compute complex shearlet coefficients (also show spectra)
[ST,Psi] = shearletTransformSpect(X,[],0);

idx = 14;

figure

subplot(3,2,1)
imagesc(X);
axis image off
colorbar
title('original image')

subplot(3,2,3)
imagesc(real(ST(:,:,idx)))
axis image off
colorbar
title('shearlet coefficients (real part)')

subplot(3,2,4)
imagesc(imag(ST(:,:,idx)))
axis image off
colorbar
title('shearlet coefficients (imaginary part)')

subplot(3,2,5)
imagesc(abs(ST(:,:,idx)))
axis image off
colorbar
title('shearlet coefficients (absolute value)')

subplot(3,2,6)
imagesc(angle(ST(:,:,idx)))
axis image off
colorbar
title('shearlet coefficients (phase)')

figure

subplot(1,3,1)
imagesc(Psi(:,:,idx))
axis image off
colorbar
title('shearlet in Fourier domain')

subplot(1,3,2)
imagesc(real(fftshift(ifft2(ifftshift(Psi(:,:,idx))))))
axis image off
colorbar
title('sherlet in time domain (real part)')

subplot(1,3,3)
imagesc(imag(fftshift(ifft2(ifftshift(Psi(:,:,idx))))))
axis image off
colorbar
title('sherlet in time domain (real part)')

%% show difference for real real shearlets
[ST,Psi] = shearletTransformSpect(X,'realReal',1);
d_real = squeeze(sum(sum(abs(imag(ST))))/256/256);

[ST,Psi] = shearletTransformSpect(X,'realReal',0);
d = squeeze(sum(sum(abs(imag(ST))))/256/256);

figure
subplot(1,2,1)
plot(d_real)
subplot(1,2,2)
plot(d)
%Note: To really see the difference you have to comment the last line in
%shearletTransformSpect.m

%% use own auxiliary function (C^\infty)
auxaux = @(x) exp(-(1./((1+x).^2) + 1./((1-x).^2)));
aux = @(x) auxaux(x-1)./(auxaux(x-1)+auxaux(x)) .*(x>=0).*(x<=1) + (x>1);

[ST,Psi] = shearletTransformSpect(X,'shearletArg',aux);

idx = 14;
figure
subplot(2,2,1)
h = imagesc(X);
axis image off
colorbar
title('original image')

subplot(2,2,2)
imagesc(ST(:,:,idx))
axis image off
colorbar
title('shearlet coefficients')

subplot(2,2,3)
imagesc(Psi(:,:,idx))
axis image off
colorbar
title('shearlet in Fourier domain')

subplot(2,2,4)
imagesc(fftshift(ifft2(ifftshift(Psi(:,:,idx)))))
axis image off
colorbar
title('sherlet in time domain')

%% use smooth shearlet
[ST,Psi] = shearletTransformSpect(X,'shearletSpect','meyerSmoothShearletSpect');

idx = 18;
figure
subplot(2,2,1)
h = imagesc(X);
axis image off
colorbar
title('original image')

subplot(2,2,2)
imagesc(ST(:,:,idx))
axis image off
colorbar
title('shearlet coefficients')

subplot(2,2,3)
imagesc(Psi(:,:,idx))
axis image off
colorbar
title('shearlet in Fourier domain')

subplot(2,2,4)
imagesc(fftshift(ifft2(ifftshift(Psi(:,:,idx)))))
axis image off
colorbar
title('sherlet in time domain')
