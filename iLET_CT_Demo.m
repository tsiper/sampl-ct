% -------------------------------------------------------------------------
% This is a demo script for the TV iLET formulation.
% The script demonstrates the TV iLET based reconstruction for the TV based problem
%
%    min_{x} 0.5||y - Hx||_2^2 + \lambda TV(x) 
%
%
% Written by Oren Solomon and Shahar Tsiper.
% -------------------------------------------------------------------------

Initialize();

% clc;
% clear;
% close all;

global DEBUG showValue updateReg VERBOSE DebugFlag;
DEBUG     = 0;
showValue = 0;
updateReg = 1;
VERBOSE   = 0;

SaveFlag = 0;
DebugFlag = 0;

%% Initialization

SNR_db = 36;                       % blurred signal to noise ratio (in dB)
Sigma = 10^-(SNR_db/20);
%% Distort (blur and add noise) original image

% The image size
N = 128;
DecFactor = 4;

x0 = LoadPhantom(N,'zubal');

% theta = (0:1/N:1-1/N)*180;
% H = paralleltomo(N,theta,N,N);

% For PSNR comparison only


% Solver parameters
% Params.L           = 0.3573;
Params.L           = 25;
% Params.lambda      = 0.001; 0.0001;
Params.lambda      = 0.02;
Params.IterMax     = 15;
Params.IterMaxIRLS = 8;

Fista_Iters = 25;

%% Defining the operators
D = DecOperator(DecFactor,'uniform');
H  = @(x) real(invF1(D(M(App(x)))));
% H  = @(x) real(invF1(D(App(x))));
Ht = @(y) real(App_T(D(F1(y))))*DecFactor;

%% Creating the measurements with additional noise

y0     = real(invF1(D(App(x0))));
y_n    = D(GaussianNoise(y0,Sigma));
psnr_0 = psnr(y_n, y0);
b = real(invF1(M(F1(y_n))));


%% Solve TV iLET
t0 = tic;

% b = real(invF1((F1(y_n))));
[ x_iLET0, obj_stack, a_stack ] = TV_iLET_CT( b, H,Ht, Params,x0 );
toc(t0);

% Calculate PSNR
psnr_1 = psnr(x_iLET0, x0);


%% Function value
figure;
% plot(obj_stack, 'linewidth', 2);xlabel('Iteration number');title('Objective function value');grid on;
x_fbp = Ht(b);
psnr_fbp = psnr(x_fbp,x0);
% x_fbp_tv = denoise_bound(x_fbp,lambda/100,0,1);

lambda_opt = fminsearch(@(lambda) -psnr(denoise_bound(x_fbp,lambda,0,1),x0),1);

x_fbp_tv0 = denoise_bound(x_fbp,lambda_opt,0,1);
x_fbp_fixed0 = denoise_bound(x_fbp,Params.lambda,0,1);

psnr_fbp_tv = psnr(x_fbp_tv0,x0);

figure;subplot(131);imagesc(x_fbp);subplot(132);imagesc(x_fbp_tv0);subplot(133);imagesc(x_iLET0);

%% Comparing to state of the art TV deblurring
% Preconditioned
% A       = @(x) real(invF1(Mnew(F1(reshape(H*vec(x),[N,M])))));
% y_Fista = real(invF1(Mnew(F1(reshape(y_n,[N,M])))));
% 
% % Non Pre-Conditioned
% A       = @(x) reshape(H*vec(x),[N,M]);
% y_Fista = reshape(y_n,[N,M]);

% At      = @(y) vec2im(H'*vec(y));


TV_Iters    = 25;
% lambda = 3;
% L = 20000;
% DebugFlag = 0;
x_fista0 = FISTA_TV_adp(b,H,Ht,zeros(N),Params.lambda,Params.L,Fista_Iters,TV_Iters,x0);
x_fista1 = FISTA_TV(b,H,Ht,zeros(N),Params.lambda,Params.L,Fista_Iters,TV_Iters,x0);

%% Comparing Images
x_iLET   = EqualizeImage(x_iLET0,x0);
x_fista  = EqualizeImage(x_fista0,x0);
x_fista_fixed  = EqualizeImage(x_fista1,x0);
x_fbp_tv = EqualizeImage(x_fbp_tv0,x0);
x_fbp_fixed = EqualizeImage(x_fbp_fixed0,x0);

%%
CompareImages(x_fista,x0,'FISTA_TV');
CompareImages(x_iLET,x0,'iLET TV');
CompareImages(x_fbp_tv,x0,'FBP Best TV');

%% Coefficients
h2=figure;
plot(1:size(a_stack, 2), a_stack, '-*');grid on;title('iLET coefficients');
semilogy(1:size(a_stack, 2), abs(a_stack), '-*');grid on;title('iLET coefficients');
xlabel('Iteration number');
% set();
legend({'x^{(n-1)}', 'x^{(n)}', '\nabla J(x^{(n)})_{TV}', 'Preconditioned \mu=1/ \tau', 'Preconditioned \mu=10/ \tau'},...
    'interpreter','tex');%, 'TV', 'Preconditioned TV \mu=1/ \tau', 'Preconditioned TV \mu=10/ \tau');

%% Plotting some nice figures

% CT Comparison
Width = 900; Height = 600;
figure('Position',[50 50 Width Height]);
g = 0.01; h = 0; w = 0;

set(0,'defaulttextinterpreter','none');
ax(1) = subtightplot(2,3,1,g,h,w);imshowzoom(x0,'GD');
ax(2) = subtightplot(2,3,2,g,h,w);imshowzoom(x_fbp_tv,'FBP');
ax(3) = subtightplot(2,3,3,g,h,w);imshowzoom(x_fbp_fixed,'FBP_Fix');
ax(4) = subtightplot(2,3,4,g,h,w);imshowzoom(x_fista_fixed,'Fista_Fix');
ax(5) = subtightplot(2,3,5,g,h,w);imshowzoom(x_fista,'Fista');
ax(6) = subtightplot(2,3,6,g,h,w);imshowzoom(x_iLET,'iLET');
linkaxes(ax);
set(0,'defaulttextinterpreter','latex');

% SaveFigure(['Reconstructions_Denoised_D',num2str(DecFactor),'_',num2str(N),...
%     '_',num2str(DesiredSNR(n)),'_F',num2str(IterFactor)],Width,Height);


%% Plotting the table
