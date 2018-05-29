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

SNR_db = 34;                       % blurred signal to noise ratio (in dB)
Sigma = 10^-(SNR_db/20);

% The image size
N = 256;
DecFactor = 3;

%% Distort (blur and add noise) original image

% PhantomName = 'zubal';
% PhantomName = 'thorax';
PhantomName = 'brain';

x0 = LoadPhantom(N,PhantomName);

% theta = (0:1/N:1-1/N)*180;
% H = paralleltomo(N,theta,N,N);

% For PSNR comparison only


% Solver parameters
% Params.L           = 0.3573;
% Params.L           = 25;
% Params.lambda      = 0.001; 0.0001;
Params.lambda      = 0.02;
Params.IterMax     = 25;
Params.IterMaxIRLS = 12;

Fista_Iters = Params.IterMax;
TV_Iters    = 100;


%% Defining the operators
D = DecOperator(DecFactor,'uniform');
H  = @(x) D(Rpp(x));
% H  = @(x) real(invF1(D(App(x))));
Ht = @(y) Rpp_T(D(M(y)))*DecFactor;

%% Finding the Lipshitz constant according to the PowerMethod
Params.L     = PowerMethod(@(x) Ht(H(x)),N);
% Params.L = 11.5; % N=256, DecFactor= 6

%% Creating the measurements with additional noise

y0     = D(Rpp(x0));
y_n    = D(GaussianNoise(y0,Sigma));
psnr_0 = psnr(y_n, y0);
b      = (y_n);

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

lambda_opt = fminsearch(@(lambda) -psnr(FGP(x_fbp,lambda,TV_Iters),x0),1);

x_fbp_tv0 = FGP(x_fbp,lambda_opt,TV_Iters);
x_fbp_fixed0 = FGP(x_fbp,Params.lambda,TV_Iters);

psnr_fbp_tv = psnr(x_fbp_tv0,x0);

%% Comparing to state of the art TV deblurring
% Preconditioned
% A       = @(x) real(invF1(Mnew(F1(reshape(H*vec(x),[N,M])))));
% y_Fista = real(invF1(Mnew(F1(reshape(y_n,[N,M])))));
% 
% % Non Pre-Conditioned
% A       = @(x) reshape(H*vec(x),[N,M]);
% y_Fista = reshape(y_n,[N,M]);

% At      = @(y) vec2im(H'*vec(y));


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
% CompareImages(x_fista,x0,'FISTA_TV');
% CompareImages(x_iLET,x0,'iLET TV');
% CompareImages(x_fbp_tv,x0,'FBP Best TV');

%% Coefficients Graph
figure('Position',[50 50 1280 720]);

PlotStyle  = {'x-.','o-','^-','*--','d-','s--','+-'};
PlotColor  = {[.8 0 0],[0 0.5 0],[.7 0 .7],[0 0 0],[0 0 .8],[0 0.5 0.5],[0.9 0.4 0]};
LineWidth  = 2;
FontSize   = 26;
FigWidth   = 1280;
FigHeight  = 720;
PlotParams = {'LineWidth',2.5,'MarkerSize',12}; 

LegendEntries = {'x^{(n-2)}','x^{(n-1)}', 'x^{(n)}', '\nabla J(x^{(n)})_{TV} - \lambda_1',...
    '\nabla J(x^{(n)})_{TV} - \lambda_2', '\nabla J(x^{(n)})_{TV} - \lambda_3',''};

g = 0;
mv = [0.17 0.03]; mh = [0.1 0.03];

% plot(1:size(a_stack, 2), a_stack, '-*');grid on;title('iLET coefficients');
% plot(SNR(SNR_Range),yPSNR_Vals(SNR_Range,i+shift),PlotStyle{i},'Color',PlotColor{i},PlotParams{:});

subtightplot(1,1,1,g,mv,mh);
for i=1:size(a_stack,1)-1
    semilogy((1:size(a_stack, 2))', abs(a_stack(i,:))', PlotStyle{i},'Color',PlotColor{i},PlotParams{:});
    hold on;
end
grid on;%title('iLET coefficients');
xlabel('Iteration number');
ylabel('iLET Coefficient Magnitude');
xlim([0 25]);
ylim([10^-4 5]);
% set();
hl = legend(LegendEntries{1:i},'location','southeast');%, 'TV', 'Preconditioned TV \mu=1/ \tau', 'Preconditioned TV \mu=10/ \tau');
set(hl,'Interpreter','tex')

set(gca,'FontSize',FontSize);

SaveFigure('iLET_Coeffs');


%% Plotting some nice figures

% CT Comparison
Width = 600; Height = 900;
figure('Position',[50 50 Width Height]);
g = 0.01; h = 0; w = 0;

x_fbp_eq = EqualizeImage(x_fbp,x0);
% x_fbp_eq = x_fbp;

set(0,'defaulttextinterpreter','none');
ax(1) = subtightplot(3,2,1,g,h,w);imshowzoom(x0,'a');
ax(2) = subtightplot(3,2,2,g,h,w);imshowzoom(x_fbp_eq,'b');
ax(3) = subtightplot(3,2,3,g,h,w);imshowzoom(x_fbp_fixed,'c');
ax(4) = subtightplot(3,2,4,g,h,w);imshowzoom(x_fista_fixed,'d');
ax(5) = subtightplot(3,2,5,g,h,w);imshowzoom(x_fista,'e');
ax(6) = subtightplot(3,2,6,g,h,w);imshowzoom(x_iLET,'f');
linkaxes(ax);
set(0,'defaulttextinterpreter','latex');

SaveFigure(['iLET_',PhantomName]);

% SaveFigure(['Reconstructions_Denoised_D',num2str(DecFactor),'_',num2str(N),...
%     '_',num2str(DesiredSNR(n)),'_F',num2str(IterFactor)],Width,Height);

%% Plotting the table
