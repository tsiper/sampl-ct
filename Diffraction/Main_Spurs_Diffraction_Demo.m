clear;
close all
clc

%%
%Optical set-up constants& parameters
INT_NUM = 74;
PrjNum = 73;
leap = 1;
lambda= 632.8*10^(-9); %wavelength in [m]
SIZE= 256;
nm=1.33; %RI of the medium- water
pix_ccd=5.2*10^(-6); %pixel size of the sensor [m]
Mag=60; %total magnification of the optical set-up
%teta = linspace(0,180,73); %73 acquisitions over 180 degrees
teta = (0:1/PrjNum:1-1/PrjNum)*180;

dX=pix_ccd/Mag; %Transversial jumps [m]
dx= dX; 

loc= [pwd];
%%
%Read
load('COM.mat') %From FIRSTRUN 
col_cm = x_c; 
row_cm = y_c; 
I= zeros(SIZE, SIZE,  INT_NUM);
temp = rgb2gray(imread(strcat(loc,'\T_cell_data\T_xy_1.png')));
I(:,:,1)= temp(row_cm(1)-SIZE/2 :row_cm(1)+ SIZE/2-1 ,col_cm(1)-SIZE/2 :col_cm(1)+ SIZE/2-1); %First cut

for a= 1:leap:INT_NUM-1
temp = rgb2gray(imread(strcat(loc,'\T_cell_data\', num2str(a),'.png')));
I(:,:,a+1)= temp(row_cm(a)-SIZE/2 :row_cm(a)+ SIZE/2-1 ,col_cm(a)-SIZE/2 :col_cm(a)+ SIZE/2-1);
end

%%
%OPDs extraction

%1st stage is turning all interferograms into OPD maps
[rows,cols]= size(I(:,:,1));
N=cols/4; 
M=rows/4;
I_1= zeros(rows, cols,  size(I,3));
I_2= zeros(M, N,  size(I,3));
I_3= zeros(rows, cols,  size(I,3));
I_4= zeros(rows, cols,  size(I,3));

for k=1: size(I,3) 
 I_1(:,:,k)= fftshift(fft2(ifftshift(I(:,:,k)))); %fftshift vs ifftshift only matters for odd grids

%cc cropping
if (k==1)
    cc_location = find_max(I_1(SIZE/2+1,:,k) , cols);
end

I_2(:,:,k) = I_1((1.5*N+1):2.5*N, cc_location-rows/8:cc_location+rows/8 -1 ,k);

%%zero padding
I_3(1.5*N+1:2.5*N,1.5*N+1:2.5*N ,k)= I_2(:,:,k);

I_3(:,:,k) = fftshift(ifft2(ifftshift(I_3(:,:,k))));
 
I_4(:,:,k)= atan2 (  imag(I_3(:,:,k))  , real(I_3(:,:,k))  );

end

N=SIZE;
M=SIZE;

A= zeros(M, N, size(I,3)-1);
%An array where the 1st cell contains the interferogram for teta= 0deg,
%the 2rd for teta- delta deg and so on. 

for l=1: size(A,3) 
A(:,:,l)= I_4(:,:,l+1) - I_4(:,:,1);

A(:,:,l)=phase_unwrapping_2D(A(:,:,l));

end

A(A<=0)=0;

%% TOMOGRPAHIC RECONSTRUCTION

%% reconstruct with radon- assume ballistics
% %reconstruct for each slice seperatly- many 2D IDFTs
 n_radon_gili= tomo_radon_2Ds(A ,lambda, pix_ccd, Mag, teta, nm );
 save('n_radon_gili','n_radon_gili');
 [n_radon_AT_SPURS,n_radon_AT_FBP] = tomo_radon_2Ds_AT(A ,lambda, pix_ccd, Mag, teta, nm );
 save('n_radon_AT','n_radon_AT_SPURS','n_radon_AT_FBP');

%reconstruct for all slices together- 3D IDFT (same result, much faster)
%n_radon= tomo_radon_3D(A ,lambda, pix_ccd, Mag, teta, nm );

%% reconstruct with diffraction- rytov
%[n_ODT, n_ODT_iter] = ODT_AT( A, I_3, lambda, pix_ccd, Mag, teta, nm);

%% Plotting AT vs. Gili
plot_n_radon(permute(n_radon_gili,[3 2 1]),permute(n_radon_AT,[3 2 1]) ,[1.33 1.4],N);

%% presentation
n_radon= (permute(n_radon, [3 2 1]));
range= [1.33 1.4];

figure
subplot(3,3,1)
imagesc(squeeze(n_radon(:,SIZE/2+1,:)),range)
title 'sagittal slice- radon'
colorbar
hold on
subplot(3,3,2)
imagesc(squeeze(n_radon(SIZE/2+1,:,:)),range)
title 'axial slice- radon'
colorbar
hold on
subplot(3,3,3)
imagesc(flipud(n_radon(:,:,SIZE/2+1)),range)
title 'coronal slice- radon'
colorbar
hold on
subplot(3,3,4)
imagesc(rot90(squeeze(n_ODT(:,SIZE/2+1,:))',2),range)
title 'sagittal slice- ODT'
colorbar
hold on
subplot(3,3,5)
imagesc(rot90(n_ODT(:,:,SIZE/2+1)',2),range)
title 'axial slice- ODT'
colorbar
hold on
subplot(3,3,6)
imagesc(fliplr(squeeze(n_ODT(SIZE/2+1,:,:))'),range)
title 'coronal slice- ODT'
colorbar
hold on
subplot(3,3,7)
imagesc(rot90((squeeze(n_ODT_iter(:,SIZE/2+1,:)))',2),range)
title 'sagittal slice- ODT with iterative'
colorbar
hold on
subplot(3,3,8)
imagesc(rot90(n_ODT_iter(:,:,SIZE/2+1)',2),range)
title 'axial slice- ODT with iterative'
colorbar
hold on
subplot(3,3,9)
imagesc(fliplr(squeeze(n_ODT_iter(SIZE/2+1,:,:))'),range)
title 'coronal slice- ODT with iterative'
colorbar