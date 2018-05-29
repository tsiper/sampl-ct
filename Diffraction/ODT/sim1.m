%Run simulation
clear all
close all
clc

%System parameters
grid_size=pow2(8);
size_x= grid_size;
size_y= grid_size; 
size_z= grid_size;
nm=1.33; %RI of the medium
ns=[1.36 1.39 1.4 1.38]; %RI of the object. Mor got 1.36 and not 1.39 for the cytplasm, but it is possibly an underestimation due to poor phase unwrapping
r=5*10^(-6);%Outer radius of the object [m]
lambda=532*10^(-9); %wavelength in [m]
delta= 5; %phase incrments for TPM
ang_range=180;
teta= 0:delta:ang_range; %control angular coverage

% Optical set-up parameters
pix_ccd=5.2*10^(-6); %pixel size of the sensor [m]
Mag=60; %total magnification of the optical set-up

dX=pix_ccd/Mag; %Transversial jumps [m]
alpha = asind((3/4)*lambda/(2*pix_ccd)); %off axis angle corrected
p_x=single(-size_x/2:(size_x/2-1));
p_y=single(-size_y/2:(size_y/2-1));
p_z=single(-size_z/2:(size_z/2-1));
[X,  Y, Z]=meshgrid(p_x*dX,p_y*dX, p_z*dX);

%create model
n_1=(((Y-0*dX)/(0.65*r)).^2+((X-45*dX)/(1*r)).^2+((Z-0*dX)/(0.45*r)).^2);
n_2=(((Y-0*dX)/(0.75*r)).^2+((X+45*dX)/(1.1*r)).^2+((Z-0*dX)/(0.45*r)).^2); 
n2=((Y+12*dX)/(0.35*r)).^2+((X+45*dX)/(0.5*r)).^2+((Z-8*dX)/(0.25*r)).^2;
n3=((Y-5*dX)/(0.3*r)).^2+((X-55*dX)/(0.45*r)).^2+((Z+4*dX)/(0.3*r)).^2;
n4=((Y-20*dX)/(0.3*r)).^2+((X+32*dX)/(0.25*r)).^2+((Z+9*dX)/(0.25*r)).^2;
n= nm*ones(size(n_1));
n(n_1<1)=ns(1); 
n(n_2<1)=ns(1); 
n(n2<=1)=ns(2);
n(n3<=1)=ns(3);
n(n4<=1)=ns(4);
clear X Y Z

% % show model slices
%  figure
%  imshow (n(:,:,grid_size/2+1),[])
%  figure
%  imshow (squeeze(n(:,grid_size/2+1,:)),[])
%  figure
%  imshow (squeeze(n(grid_size/2-1,:,:)),[])

% create interfogram array &axtract the phase from it. 1st is with no sample
% sample rotation
I = INTERFEROGRAM_GENERATOR_imrotate(n,teta, lambda, alpha, nm, dX, pix_ccd);

%%
%allocate space
SIZE=grid_size;
[rows,cols,num]= size(I); 
N=cols/4; 
M=rows/4;
I_1= zeros(rows, cols,  num);
I_2= zeros(M, N,  num);
I_3= zeros(rows, cols,  num);
I_4= zeros(rows, cols,  num);
A= zeros(rows, cols,  num-1);
Miguel= zeros(rows, cols,  num-1);

%% extract phase with algorithm A
for k=1: num 
 I_1(:,:,k)= fftshift(fft2(ifftshift(I(:,:,k)))); %fftshift vs ifftshift only matters for odd grids

max_1=0;
%cc cropping
if (k==1) %off axis angle remains constant
    for l=1:SIZE
    [cc_col, max_row] = find_max2(I_1(l,:,k) , cols);
    if max_row>max_1
        max_1= max_row;
        cc_location_row =l; 
        cc_location_col= cc_col;
    end
    end
end
  
I_2(:,:,k) = I_1(cc_location_row-SIZE/8:cc_location_row+SIZE/8 -1, cc_location_col-SIZE/8:cc_location_col+SIZE/8 -1 ,k);

%%zero padding
I_3(1.5*N+1:2.5*N,1.5*N+1:2.5*N ,k)= I_2(:,:,k);

I_3(:,:,k) = fftshift(ifft2(ifftshift(I_3(:,:,k)))); %complex field
 
I_4(:,:,k)= atan2 (  imag(I_3(:,:,k))  , real(I_3(:,:,k))  );

end

for k=1: num-1   
A(:,:,k)= Miguel_2D_unwrapper(single(I_4(:,:,k+1) - I_4(:,:,1))); %phase
end

% %% alternatively: calculate the phase straight from the model
% L= length(teta);
% n_dif= n-nm;
% 
% phase_true= zeros(size_x,size_y,L);
% for k=1:L
% for t=1:size_x 
% rot= imrotate((squeeze(n_dif(t,:,:))'), teta(k) , 'nearest','crop'); 
% phase_true(t,:,k)= (sum(rot,1)*dX*2*pi/lambda);
% end
% end


%% TOMOGRPAHIC RECONSTRUCTION

%% reconstruct with radon- assume ballistics
PolarGrid = BuildPolarGrid(256,73);
% %reconstruct for each slice seperatly- many 2D IDFTs
n_radon= tomo_radon_2D(A ,lambda, dX, teta, nm, PolarGrid);

%reconstruct for all slices together- 3D IDFT (same result, much faster)
%n_radon= tomo_radon_3D(A ,lambda, pix_ccd, Mag, teta, nm );

%% reconstruct with diffraction- rytov
[n_ODT, n_ODT_iter] = ODT( A, I_3, lambda, pix_ccd, Mag, teta, nm);

%% presentation
n= (permute(n, [3 2 1]));
n_radon= (permute(n_radon, [3 2 1]));
range= [1.33 1.4];

figure
subplot(4,3,1)
imagesc(squeeze(n_radon(:,SIZE/2+1,:)),range)
title 'sagittal slice- radon'
colorbar
hold on
subplot(4,3,2)
imagesc(squeeze(n_radon(SIZE/2+1,:,:)),range)
title 'axial slice- radon'
colorbar
hold on
subplot(4,3,3)
imagesc(flipud(n_radon(:,:,SIZE/2+1)),range)
title 'coronal slice- radon'
colorbar
hold on
subplot(4,3,4)
imagesc(rot90(squeeze(n_ODT(:,SIZE/2+1,:))',2),range)
title 'sagittal slice- ODT'
colorbar
hold on
subplot(4,3,5)
imagesc(rot90(n_ODT(:,:,SIZE/2+1)',2),range)
title 'axial slice- ODT'
colorbar
hold on
subplot(4,3,6)
imagesc(fliplr(squeeze(n_ODT(SIZE/2+1,:,:))'),range)
title 'coronal slice- ODT'
colorbar
hold on
subplot(4,3,7)
imagesc(rot90((squeeze(n_ODT_iter(:,SIZE/2+1,:)))',2),range)
title 'sagittal slice- ODT with iterative'
colorbar
hold on
subplot(4,3,8)
imagesc(rot90(n_ODT_iter(:,:,SIZE/2+1)',2),range)
title 'axial slice- ODT with iterative'
colorbar
hold on
subplot(4,3,9)
imagesc(fliplr(squeeze(n_ODT_iter(SIZE/2+1,:,:))'),range)
title 'coronal slice- ODT with iterative'
colorbar
hold on
subplot(4,3,10)
imagesc((squeeze(n(:,SIZE/2+1,:))),range)
title 'sagittal slice- true RI'
colorbar
hold on
subplot(4,3,11)
imagesc(squeeze(n(SIZE/2+1,:,:)),range)
title 'axial slice- true RI'
colorbar
hold on
subplot(4,3,12)
imagesc(n(:,:,SIZE/2+1),range)
title 'coronal slice- true RI'
colorbar
