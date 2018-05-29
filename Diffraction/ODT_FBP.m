clear;
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

%%
t = 128;
Sinogram = squeeze(A(128,:,:));

figure;
subplot(121);imagesc(real(Sinogram));
subplot(122);imagesc(imag(Sinogram));

%% Filtered back projection

x_fbp = iradon(Sinogram,teta,N);
% x_fbp = imresize(x_fbp,[N N]);
x_fbp = x_fbp / max(x_fbp(:));

figure;
imagesc(x_fbp(:,:,end));title('FBP');
imagesc(imrotate(x_fbp(:,:,end),180));title('FBP - 180');

% %% Running with SPURS
% PolarGrid = BuildPolarGrid(size(Sinogram,1),teta);
% PhantomSamples = -vec(F1(Sinogram));
% 
% SPURS_Params = SPURS_DefaultParams(N);
% 
% OutputImages = SPURS(PhantomSamples,PolarGrid,SPURS_Params);
% 
% x_spurs = OutputImages(:,:,end);
% x_spurs = x_spurs / max(x_spurs(:));
% 
% %% Plotting SPURS
% figure;
% subplot(131);imagesc(OutputImages(:,:,end));title('SPURS');
% subplot(132);imagesc(x_fbp(:,:,end));title('FBP');
% CompareImages(x_spurs,ShiftImage(x_fbp,x_spurs));
% colormap('bone');