%% Initializing
 
clear;
Initialize;
close all;
clc;
 
%% Ct Specifications
M=177; %projections number, maximum - 179 , mod(M,down_freq) = 0;
num_of_rows=64; %in the detectors array
detectors_in_row=750; %typical number move between 750 and 800;
detector_array_y_size=1.25*10^-3; %meters
detector_array_x_size=detector_array_y_size; %an assumption;
omega=2*pi*(10/3); %an assumption of the angular velocity of the machine;
v_b=64*10^-3; %m/sec, typical speed of the table can be 26.67,48 or 64*10^-3;
R=85*10^-2; %the internal radius in meters, the outside is about 150*10^-2;
theta=floor((0:1/(M-1):1)*179);
 
%% Loading the phantom
 
N = 256;
x0=LoadPhantom(N,'brain');
imagesc(x0);
 
%% functions of boddy (x,y,z) point movement - not in use
 
x=@(t) R*cos(omega*t);
y=@(t) R*sin(omega*t);
z=@(t) v_b*t;
figure(1)
h1=fplot3(x,y,z);
title('3D description of (x,y,z) boddy point on the detectors array');
figure(2)
h2=fplot(x,z);
title('2D description of (x,y,z) boddy point on the detectors array');
figure(3)
h3=fplot(y,z);
title('2D description of (x,y,z) boddy point on the detectors array');
 
%% Designing a blurring matrix
Blurring_Level=0;
D=0.5*diag(ones(length(x0),1));
for B=1:1:Blurring_Level
    D = D+ diag(ones(length(x0)-B,1),B);
end
D = D+D';
figure(4)
subplot(321);imagesc(x0);
subplot(322);imagesc(D*x0);
subplot(323);imagesc(x0*D);
subplot(324);imagesc(D*x0*D);
x_blur=D*x0*D;
 
%% Insert Gaussian Noise
Noise_Level=0;
x_blur_noise=awgn(x_blur,1/Noise_Level,'measured');
subplot(325);imagesc(x_blur_noise);
 
%% Sinogram and Radon Transform
 
figure(5)
y_radon_original=radon(x0,theta);
y_radon_blur_noise=radon(x_blur_noise,theta);
subplot(121);imagesc(y_radon_original);
subplot(122);imagesc(y_radon_blur_noise);
[row_y_radon,col_y_radon]=size(y_radon_blur_noise);
 
%% simulating shake
%DecFactor = R*omega*detector_array_y_size/(v_b*detector_array_x_size);
%DecFactor = ceil(DecFactor/detectors_in_row*N);
DecFactor=3;
DecMethod = 'detector';
down_freq=round(M/num_of_rows);
 
D = DecOperator(DecFactor,DecMethod,down_freq);
 
y_hat = D(y_radon_blur_noise);
 
%D = DecOperator(DecFactor,DecMethod,1);
 
%Dup - Duplicates One LR matrix columns so it's size will be like HR matrix's size 
Dup=zeros([M/down_freq,M]);
size_Dup=size(Dup);
counter=0;
for i=1:size_Dup(1)
   for j=1:size_Dup(2)
        if (j > counter && j<=counter + down_freq)
            Dup(i,j) = 1;
        end
    end
    counter = counter + down_freq;
end
 
All_W = sparse([]);
All_LR_col=[];
for i=1:down_freq
    %Low - reduce HR matrix columns so it's size will be like LR matrix's size 
    Low=zeros([M,M/down_freq]);
    size_Low=size(Low);
    counter=0;
    for col=1:size_Low(2)
        for row=1:size_Low(1)
            if (row == counter*down_freq + i)
                Low(row,col) = 1;
            end
        end
        counter = counter + 1;
    end
    F=@(y) D(y)*Low*Dup;
    One_LR_col=y_hat(:,i:down_freq:M)*Dup;
    All_W=[All_W;BuildWfromD(F,size(y_hat))];
    One_LR_col = One_LR_col(:);
    All_LR_col = [All_LR_col ; One_LR_col];
end

A = @(x) All_W*x;
At =@(x) All_W'*x;
%% running FISTA
 
x_start=vec(zeros(size(y_radon_original)));
x_out= FISTA_L2( All_LR_col, A, At, x_start, 0,20, 300);
%x_out= FISTA_TV( All_LR_col,A,At,x_start,1e-3,100,1000,20,y_radon_original);
x_out_recon=reshape(x_out,size(y_radon_original));
%% check PSNR
 
PSNR_sinogram = psnr(x_out_recon , y_radon_original);
figure(6)
subplot(121);imagesc(y_radon_original);
subplot(122);imagesc(x_out_recon);
PSNR_phantom = psnr(iradon(y_radon_original,theta) , iradon(x_out_recon,theta));
figure(7)
subplot(121);imagesc(iradon(y_radon_original,theta));
subplot(122);imagesc(iradon(x_out_recon,theta));