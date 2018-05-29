clc;
clear;
close all;

%% Generate image - Denoising
% -------------------------------------------------------------------- %
X = double(imread('cameraman.pgm'));
X = X/255;

% Add noise
Nsig = 2e-2;
randn('seed',314);
Bobs = X + Nsig*randn(size(X));


% Run solution
% -------------------------------------------------------------------- %
lambda = 0.02;

% Denoising
X_den=denoise_bound(Bobs,lambda,0,Inf);

% Show images
figure;
h(1) = subplot(2,3,1);
imshow(X,[]);title('Ground truth');
h(2) = subplot(2,3,2);
imshow(Bobs,[]);title(['Noisy, \sigma = ' num2str(Nsig)]);
h(3) = subplot(2,3,3);
imshow(X_den,[]);title('TV denoised');


%% Generate image - Deblurring
% -------------------------------------------------------------------- %
X = double(imread('cameraman.pgm'));
X = X/255;


P=1/9*ones(7,7);
center=[4,4];

Nsig = 1e-4;
randn('seed',314);
Bobs = imfilter(X,P,'symmetric')+Nsig*randn(size(X));

% Run solution
% -------------------------------------------------------------------- %
lambda = 1e-3;
pars.fig = 0;

X_deblur=deblur_tv_fista(Bobs,P,center,lambda,0,Inf, pars);

% Show images
h(4) = subplot(2,3,4);
imshow(X,[]);title('Ground truth');
h(5) = subplot(2,3,5);
imshow(Bobs,[]);title(['Blurred and noisy, \sigma = ' num2str(Nsig)]);
h(6) = subplot(2,3,6);
imshow(X_deblur,[]);title('TV deblurred');
linkaxes(h, 'xy');

%% Blur
% -------------------------------------------------------------------- %
M = 128;
Band = 10;
BlurWidth = 2;0.8;

%% Generate blurred image
[A,b,x] = blur(M, Band ,BlurWidth);
Bobs = reshape(b, M, M);

% Add noise
Nsig = 2e-2;
randn('seed',314);
Bobs = Bobs + Nsig*randn(size(Bobs));

psf = reshape(full(A(:, 10000)), 128, 128);
psf_f = fft2(psf);
psf_fc = abs(psf_f);
psf_c = ifftshift(ifft2(psf_fc));
psf_center = [65 65];

% Run solution
% -------------------------------------------------------------------- %
lambda = 1e-4;
pars.fig = 0;

X_deblur=deblur_tv_fista(Bobs,psf_c,psf_center,lambda,0,Inf, pars);

figure;
hb(1) = subplot(2,3,1);
imshow(reshape(x, M, M),[]);title('Ground truth');
hb(2) = subplot(2,3,2);
imshow(Bobs,[]);title(['Blurred and noisy, \sigma = ' num2str(Nsig)]);
hb(3) = subplot(2,3,3);
imshow(X_deblur,[]);title('TV deblurred');
linkaxes(hb, 'xy');
subplot(2,3,4);
imshow(psf_c, []);title('PSF');
















