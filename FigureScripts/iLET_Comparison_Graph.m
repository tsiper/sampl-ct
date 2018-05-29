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

SaveFlag = 0;
DebugFlag = 0;

%% Initialization

SNR_db = 20:50;                       % blurred signal to noise ratio (in dB)
Sigma = 10.^-(SNR_db/20);

% The image size
N = 128;
DecFactor = 3;

%% Distort (blur and add noise) original image
x0 = LoadPhantom(N,'zubal');

% For PSNR comparison only

% Solver parameters
Params.lambda      = 0.02;
Params.IterMax     = 25;
Params.IterMaxIRLS = 12;

Fista_Iters = Params.IterMax;
TV_Iters    = 100;

%% Defining the operators
D  = DecOperator(DecFactor,'uniform');
H  = @(x) D(Rpp(x));
Ht = @(y) Rpp_T(D(M(y)))*DecFactor;

%% Finding the Lipshitz constant according to the PowerMethod
Params.L = PowerMethod(@(x) Ht(H(x)),N);
% Params.L = 11.5; % N=256, DecFactor= 6

%% Creating the initial measurements with additional noise
y0     = D(Rpp(x0));

% Opening the waitbar
wb = MyWaitbar(0,'Computing all SNRs for iLET Graph');

% Iterating for various noise levels
for n = 1:length(Sigma)
    % Adding noise to the measurements
    y_n    = GaussianNoise(D(y0),Sigma(n));

    % Running a single iteration, equivalent to filtered back projection
    x_fbp = Ht(y_n);
    
    % Finding the optimal lambda by exhaustive grid search
    lambda_opt   = fminsearch(@(lambda) -psnr(FGP(x_fbp,lambda,TV_Iters),x0),1);
    x_fbp_tv0    = FGP(x_fbp,lambda_opt,TV_Iters);
    x_fbp_fixed0 = FGP(x_fbp,Params.lambda,TV_Iters);

    % Solve TV iLET
    [ x_iLET0, obj_stack, a_stack ] = TV_iLET_CT( y_n, H,Ht, Params,x0 );

    % Solving with FISTA for comparison
    x_fista0 = FISTA_TV_adp(y_n,H,Ht,zeros(N),Params.lambda,Params.L,Fista_Iters,TV_Iters,x0);
    x_fista1 = FISTA_TV(y_n,H,Ht,zeros(N),Params.lambda,Params.L,Fista_Iters,TV_Iters,x0);

    %% Equalizing all the images for fair comparison
    x_iLET         = EqualizeImage(x_iLET0,x0);
    x_fista        = EqualizeImage(x_fista0,x0);
    x_fista_fixed  = EqualizeImage(x_fista1,x0);
    x_fbp_tv       = EqualizeImage(x_fbp_tv0,x0);
    x_fbp_fixed    = EqualizeImage(x_fbp_fixed0,x0);

    %% PSNR Calculations
    PSNR_fbp_opt(n)   = psnr(x_fbp_tv,x0); %#ok<*SAGROW>
    SSIM_fbp_opt(n)   = ssim(x_fbp_tv,x0);
    
    PSNR_fbp_fix(n)   = psnr(x_fbp_fixed,x0);
    SSIM_fbp_fix(n)   = ssim(x_fbp_fixed,x0);
    
    PSNR_fista_fix(n) = psnr(x_fista_fixed,x0);
    SSIM_fista_fix(n) = ssim(x_fista_fixed,x0);
    
    PSNR_fista_opt(n) = psnr(x_fista,x0);
    SSIM_fista_opt(n) = ssim(x_fista,x0);
    
    PSNR_ilet(n)      = psnr(x_iLET,x0);
    SSIM_ilet(n)      = ssim(x_iLET,x0);
    
    PSNR_Sinogram(n) = psnr(y_n, y0);
    
    % Updating waitbar
    MyWaitbar(n/length(Sigma),wb);
    

    
end
close(wb); % closing the waitbar


%% Setting up plot variables
LegendEntries = {...
    'FBP+TV - Constant $\lambda$',...
    'FBP+TV - Optimal $\lambda$',...
    'FISTA+TV - Constant $\lambda$',...
    'FISTA+TV - Adaptive $\lambda$',...
    'TV-iLET  - Adaptive $\lambda$','location','southeast'};
PSNR_cell = {PSNR_fbp_fix,PSNR_fbp_opt,PSNR_fista_fix,PSNR_fista_opt,PSNR_ilet};
SSIM_cell = {SSIM_fbp_fix,SSIM_fbp_opt,SSIM_fista_fix,SSIM_fista_opt,SSIM_ilet};

%% Plot specifications
PlotStyle  = {'x-.','o-','^-','*--','d-','s--'};
PlotColor  = {[.8 0 0],[0 0.5 0],[.7 0 .7],[0 0 0],[0 0 .8],[0 0.5 0.5]};
LineWidth  = 2;
FontSize   = 32;
FigWidth   = 1280;
FigHeight  = 720;
PlotParams = {'LineWidth',2.5,'MarkerSize',12}; 

g = 0;
mv = [0.17 0.03]; mh = [0.1 0.03];

%% Plotting

% Latex Interpreter

% PSNR Graph
figure('Position',[50 50 FigWidth FigHeight]);
subtightplot(1,1,1,g,mv,mh);
for i=1:5
    plot(PSNR_Sinogram,PSNR_cell{i},PlotStyle{i},'Color',PlotColor{i},PlotParams{:}); hold on;
end
legend(LegendEntries{:}); grid on;
set(0,'defaulttextinterpreter','latex');
xlabel('Input Sinogram SNR'); ylabel('Reconstruction PSNR [dB]');
set(gca,'fontsize',FontSize);

SaveFigure('iLET_PSNR_Graph');

% SSIM Graph
figure('Position',[50 50 FigWidth FigHeight]);
subtightplot(1,1,1,g,mv,mh);
for i=1:5
    plot(PSNR_Sinogram,SSIM_cell{i},PlotStyle{i},'Color',PlotColor{i},PlotParams{:}); hold on;
end
legend(LegendEntries{:}); grid on;
xlabel('Input Sinogram SNR [dB]'); ylabel('Reconstruction SSIM');
set(gca,'fontsize',FontSize);

SaveFigure('iLET_SSIM_Graph');