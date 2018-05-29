%% try to combine cp_als in functional 

%% In case of emergency press ctrl+enter
clear;clc;close all;
%% Global parameters and Files
global DebugFlag;
DebugFlag = 0;
% when working with the files which are contained in this folder we use
% addpath(genpath(pwd));
% when working with shahar's files we include the fileparts inorder to gain the parent folder directory
addpath(genpath(fileparts(pwd)));

%% Creating Sampled Object
original_object = get_3D_im('Functional/SimulatorFiles'); 
%time sample rate 10^-1 sec - continuous. sample grid for Original_object
%is [0 0.1 0.2 ... 49] (491 slices).

T=10;%%space between sinograms in index
num_pix = 256;

% time vector with which we are working
dt = 1; %sample rate in second
% dt = 5; %sample rate in second
time = 0:dt:49;
index_time = int16(T*time+1);

sampled_object = original_object(:,:,index_time);
sampled_object_non_normalized=sampled_object;

%normaliation for mfista- doesnt work without it
sampled_object=sampled_object/max(sampled_object(:));
sampled_object_2d = reshape(sampled_object,num_pix^2,length(time));
clear('original_object');


%% Creating Sinograms

A_op  = @(x) real(invF1(App(x)));
At_op = @(y) real(App_T(M(F1(y))));

D = eye(length(time)) + diag(repmat(-1,1,length(time)-1),-1);
D_op = @(x) D*x';
Dt_op = @(x) D'*x;
DtD   = @(x) Dt_op(D_op(x))';

Final_A=@(x) App_2d(x,A_op);
Final_At=@(x) App_T_2d(x,At_op,num_pix);

% Getting the AtA operator
[ AtA ] = BuildAtA( num_pix ,2,'uniform' );

%% sinograms

% Y = Aradon(sampled_object_2d);
Sigma   = 0.05; % The noise level
% Sigma   = 0.00; % The noise level


% dtheta = 0.5;                             % How many degrees we jump in each scan
dtheta = 2;
theta = -90:dtheta:180;                 % The angular range - this should be bigger than -45:135
At_Y = zeros(num_pix^2,length(time));
DecFactor = 2;
DecOp = DecOperator(DecFactor,'uniform');
for t=1:length(time)
    y   = radon(sampled_object(:,:,t),theta);               % The radon projection
    %[y_n,SNR] = GaussianNoise(y,Sigma); % Adding noise
    y_n = y + Sigma.*range(y).*randn(size(y));
    y_pp = InterpolateSinogram(y_n,2*num_pix,theta,'spline');
    
    At_y = real(App_T(DecOp(M(F1(y_pp)))));
    % Trimming unnecessary zeros
    
    At_y = At_y(num_pix/2+1:3*num_pix/2,num_pix/2+1:3*num_pix/2);
    At_Y(:,t)    = vec(At_y);
end
lip = PowerMethod(AtA,num_pix);

AtA_2d_op=@(x) AtA_2d(x,AtA);

%% from here we begin

%% creating a tensor from the At_y

% convert At_Y to a 3 dimantion double
At_Y_3dim = reshape(At_Y,num_pix,num_pix,length(time));
At_Y_tensor = tensor(At_Y_3dim);
number_of_ranks = 100;

%% Collect Data for Compare Between The Ground Truth and the Projections with several CP RANK.

% collect the data of the full image (all the projections) with many CP rank
At_Y_Tensor_MultCP = cell(1,number_of_ranks);

for rank = 1:number_of_ranks
    At_Y_Tensor_MultCP{1,rank} = double(cp_als(At_Y_tensor,rank));
end

 %% find the best rank for CP

%sampled_object_cell = cell(sampled_object);
for rank = 1:number_of_ranks
    % Collect data of the PSNR between the Tensor+Noise+CP and the clean Tensor
    PSNR_vec_ALL(rank) = psnr(At_Y_Tensor_MultCP{1,rank},sampled_object);
end
%find the max PSNR (value and index)
[maxPSNRval_ALL , CP_Rank] = max(PSNR_vec_ALL(:));    

At_Y_Tensor_BestCP = At_Y_Tensor_MultCP{CP_Rank};

At_Y_CP = reshape(At_Y_Tensor_BestCP,num_pix^2,length(time));

%% from here we need to Reconstruct the At_Y_CP, and compare to the original functional.

%% Using Fista- Acheiving Reconstruction
lambda1=4*10^-3;
lambda2 = 10^-3;

TV_iters = 25;
iter = 20; % FISTA Iteration

lip = 2*lip;

%X0 = zeros(num_pix,num_pix,length(time));
X0 = zeros(num_pix^2,length(time));
%X_hat = zeros(num_pix,num_pix,length(time));

tic
%for i=1:length(time)
    %X_hat(:,:,i)=FISTA_TV(Y(:,:,i),A_op,At_op,X0(:,:,i),lambda,lip,iter,iter);
%     FISTA_functional( b,D,Dt,DtD,A,At,AtA,x_start,lambda1,lambda2,L,iters,TV_Iters,x_true)
%     DebugFlag = 1;
    X_hat=FISTA_functional(At_Y,D_op,Dt_op,DtD,Final_A,Final_At,AtA_2d_op,X0,lambda1,lambda2,lip,iter,TV_iters);
    
    X_hat_CP=FISTA_functional(At_Y_CP,D_op,Dt_op,DtD,Final_A,Final_At,AtA_2d_op,X0,lambda1,lambda2,lip,iter,TV_iters);
%end
toc


%% Normalization and reshaping to 3D
X_hat_non_normalized=X_hat;
X_hat_non_normalized = reshape(X_hat_non_normalized,num_pix,num_pix,length(time));
%normalization- NOT SURE IF NEEDED OR NOT- I THINK NOT
X_hat=X_hat/max(X_hat(:));
X_hat = reshape(X_hat,num_pix,num_pix,length(time));

%% CP -------- Normalization and reshaping to 3D
X_hat_non_normalized_CP=X_hat_CP;
X_hat_non_normalized_CP = reshape(X_hat_non_normalized_CP,num_pix,num_pix,length(time));
%normalization- NOT SURE IF NEEDED OR NOT- I THINK NOT
X_hat_CP=X_hat_CP/max(X_hat_CP(:));
X_hat_CP = reshape(X_hat_CP,num_pix,num_pix,length(time));

%% Calculating ssim and psnr for X_hat
psnr_X_hat=psnr(X_hat,sampled_object)
% ssim_X_hat=ssim(X_hat,sampled_object)

%% CP ---------------- Calculating ssim and psnr for X_hat
psnr_X_hat_CP=psnr(X_hat_CP,sampled_object)
% ssim_X_hat_CP=ssim(X_hat_CP,sampled_object)

%% Video Comparing checkup
mat2gray_scale_sampled_object=30;
fix_scale=1;

Video_checkup( real(X_hat_non_normalized), sampled_object_non_normalized, time, mat2gray_scale_sampled_object, fix_scale )

Video_checkup( real(X_hat_non_normalized_CP), sampled_object_non_normalized, time, mat2gray_scale_sampled_object, fix_scale )

%% Deriving C
X_2d = reshape(X_hat,num_pix^2,length(time));
%changing Xhat to N^2,t

[ C_hat ] = Get_C( X_2d );

%% CP ----------- Deriving C
X_2d_CP = reshape(X_hat_CP,num_pix^2,length(time));
%changing Xhat to N^2,t

[ C_hat_CP ] = Get_C( X_2d_CP );

%% Comparing C
Compare_C(time,C_hat)
%%
Compare_C(time,C_hat_CP)

%% Definition of Class- Cheating
[ class ] = get_class();
%% Using Gauss Newton to acheive K according to h(t) model
class_reshape=class(:);

%parameters
precision = 1e-5; %stop condition
itermax = 100; %number of iterations
epsilon = 10^-15; % Precision of approximated derivation
A0 = [5,1,0.1];
print=0;
epsilon_inv=10^-10;

tic
[Aopt, R, K, sum_fail, sum_fail2, sum_fail3]=Gauss_Newton_per_image( time,C_hat,dt,X_2d,A0,itermax,epsilon,precision,print,class_reshape,epsilon_inv);
[Aopt_CP, R, K_CP, sum_fail, sum_fail2, sum_fail3]=Gauss_Newton_per_image( time,C_hat_CP,dt,X_2d_CP,A0,itermax,epsilon,precision,print,class_reshape,epsilon_inv);
toc
%% Checking K
figure(100);
plot(K)
%%
figure(105);
plot(K_CP)

%% Comparing parameters MTT CBV CBF TTP and calculating ssim and psnr
[psnr_vec, ssim_vec]=Compare_Parameters( Aopt, X_2d, class, dt)
%%
[psnr_vec_CP, ssim_vec_CP]=Compare_Parameters( Aopt_CP, X_2d_CP, class, dt)

%% EXTRA WORK FOR SANITY CHECK ONLY

% % Reconstructing X2
% [ X2 ] = Deriving_X_from_K( C_hat, K, time );
% 
% % Calculating ssim and psnr for X2
% psnr_X2=psnr(X2,sampled_object)
% ssim_X2=ssim(X2,sampled_object)
% 
% % video checkup for X2- SANITY CHECK
% mat2gray_scale_sampled_object=30;
% fix_scale=1;
% 
% Video_checkup( X2, sampled_object_non_normalized, time, mat2gray_scale_sampled_object, fix_scale )
