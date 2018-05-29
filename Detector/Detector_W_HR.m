%% Initializing

clc
clear
close all

N = 32;
Mproj  = 2*(N+1);%DecFactor*(N+1)*2;
sigma_x = 0.0;       % Image noise
sigma_y = 0.0;       % Measurement noise
SNR = 300;
PhantomType = 'brain';
DecFactor = 4;

%% Decimation parameters
%DecMethod = 'uniform';
% DecMethod = 'random';
% DecMethod = 'limited';
%DecMethod = 'detector';
%DecMethod = 'pattern algorithm';
DecMethod = 'detector_mult_theta';

x0 = LoadPhantom(N,PhantomType);
Wx = 1;
Wy = 50;

D = DecOperator(DecFactor,64);

x0_noise = x0 + sigma_x * range(x0(:)) * randn(size(x0));

%% FBP and PP calculations
% s       = (-1+4/Mproj):(4/Mproj):1;
theta   = (-1/4:1/Mproj:3/4-1/Mproj)*180;
% theta_pp   = [atan(s),atan(s)+pi/2]*180/pi;
y_radon = radon(x0_noise,theta);
[row_y_radon, col_y_radon] = size(y_radon);
y_hat = real(invF1(D(F1(y_radon))));

W_1 = (y_hat(:)*y_radon(:)' + eye(row_y_radon*col_y_radon)*1e-5)/(y_radon(:)*y_radon(:)' + eye(row_y_radon*col_y_radon)*1e-5);

y_radon_2 = circshift(y_radon, [0,1]);
y_hat_2 = real(invF1(D(F1(y_radon_2))));
y_hat_2 = circshift(y_hat_2, [0,-1]);

W_2 = (y_hat_2(:)*y_radon(:)' + eye(row_y_radon*col_y_radon)*1e-5)/(y_radon(:)*y_radon(:)' + eye(row_y_radon*col_y_radon)*1e-5);

y_radon_3 = circshift(y_radon, [0,2]);
y_hat_3 = real(invF1(D(F1(y_radon_3))));
y_hat_3 = circshift(y_hat_3, [0,-2]);

W_3 = (y_hat_3(:)*y_radon(:)' + eye(row_y_radon*col_y_radon)*1e-5)/(y_radon(:)*y_radon(:)' + eye(row_y_radon*col_y_radon)*1e-5);

y_radon_4 = circshift(y_radon, [0,3]);
y_hat_4 = real(invF1(D(F1(y_radon_4))));
y_hat_4 = circshift(y_hat_4, [0,-3]);

W_4 = (y_hat_4(:)*y_radon(:)' + eye(row_y_radon*col_y_radon)*1e-5)/(y_radon(:)*y_radon(:)' + eye(row_y_radon*col_y_radon)*1e-5);

%% Running fista
 all_W=[W_1;W_2;W_3;W_4];
 all_LR_col=[y_hat(:);y_hat_2(:);y_hat_3(:);y_hat_4(:)];
 x_start=y_hat;
field1 = 'L0';  value1 = 10*max(eig(all_W'*all_W));
field2 = 'eta';  value2 = 1.1;
field3 = 'Lambda';  value3 = 0;
field4 = 'IterMax';  value4 = 2000;
field5 = 'Type';  value5 = "const";
field6 = 'NonNegOrth';  value6 = 1;
Params = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6);
[ x_out, x_history ]=FISTA(all_W,all_LR_col,zeros(1280,1280),x_start(:),Params);
x_out_recon=reshape(x_out,row_y_radon,col_y_radon);