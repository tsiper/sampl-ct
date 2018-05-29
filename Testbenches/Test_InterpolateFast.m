%% Resampling Radon to PP using our specific sub-space Kernel

%% Initializing
Initialize;
DebugFlag = 0;
tic;

%% PARAMETERS

% Chossing a phantom
PhantomTypes = {'brain','thorax','shepp','zubal'};
PhantomType  = PhantomTypes{4};

% The interpolation methods to run
% InterpMethod = {'filt'};
% InterpMethod = {'nearest','linear','pchip','spline','filt'};
% InterpMethod = {'spline','filt'};
InterpMethod = {'spline','filt'};

% Main run paramters
N          = 256;            % The image side resolution
% WindowSize =  10;            % The radius of the window
% SNR = 30;              % Signal to noise ratio
% DesiredSNR = 25:3:65;
% DesiredSNR = 30:5:35;
DesiredSNR = 80;
Sigma      = 10.^(-DesiredSNR/20); % Noise variance
% Sigma      = 0;             % Noise variance
Mp         = 2*N + 2;       % Number of projection angles
cg_iters   = 3;             % Maximum iterations for conjugate gradient alg.
DecFactor  = 1;             % Decimation factor of the measurements
% ScanType   = 'MATLAB';      % The scanning approach {'MATLAB' / 'AIR'}

% The default parameters for the kernel
B         = (Mp/10);
% R         = 1/pi/2/sqrt(2);
R         = (1/pi/5);
% R         = 1/pi/2;
W         = (pi)*(2*N+1);

% TV Fista Params
% PP Params
lambda=6e-3; % Good for 128
% lambda = 0.05;
lambda_fbp = 3;
L = 35;
FISTA_Iters = 30;
TV_Iters = 25;

%% The results matrix initialization
% ResultsRows  = InterpMethod;
ResultsRows  = [InterpMethod,cellfun(@(x) [x,'-d'],InterpMethod,'UniformOutput',0)];
ResultsRows  = [ResultsRows,{'fbp+tv','spurs+tv','art','art+tv','polar+tv'}];
ResultsCols  = {'yPSNR','ySSIM','xPSNR','xSSIM','tvPSNR','tvSSIM'};
Results      = zeros(length(ResultsRows),length(ResultsCols),length(Sigma));

%% Initializing Variables and cells

% The Pseudo-Polar singoram (y) and reconstructions (x)
x_fbp  = cell(1,length(Sigma));
x_rd_tv = x_fbp; x_art = x_fbp; x_art_tv = x_fbp; 
x_spurs = x_fbp; x_spurs_tv = x_fbp; x_fbp_tv = x_fbp;
y_pp   = cell(length(InterpMethod),length(Sigma));
x_pp   = y_pp; 
y_pp_d = y_pp; 
x_pp_d = y_pp; 
x_tv_d = y_pp;
x_tv = y_pp; 

% Loading the Phantom
x0 = LoadPhantom(N,PhantomType);

% Calculating the ideal Pseudo-Polar sinogram
y_pp0 = real(invF1(App(x0)));

% The real SNR and PSNR vectors
PSNR = zeros(size(Sigma));
SNR  = zeros(size(Sigma));

%% Generating the scanned data
% The angle ranges
N_rd     = 2*ceil(sqrt(2)*(N-floor((N-1)/2)-1))+3; % Number of sensors
theta      = (-1/4:1/Mp:3/4-1/Mp)*180;             % Radon Angles
theta_pad  = (-1/2:1/Mp:1-1/Mp)*180;               % Padded Angles for Radon
% The desired PP scanning angles
l = -N/2:N/2; theta_pp = [atan(2*l/N),fliplr(pi/2-atan(2*l/N))]*180/pi;

% The projections
y_rd     = radon(x0,theta);
y_rd_pad = radon(x0,theta_pad);
y_rd_pad = padcols(y_rd_pad,2*N+1);

% Decimating acquisition angles
theta     = theta(1:DecFactor:end);
theta_pad = theta_pad(1:DecFactor:end);

% FISTA-TV with Polar FBP Operators
A_rd = paralleltomo(N,theta,N_rd);
A    = @(x) reshape(A_rd*x(:),N_rd,length(theta));
At   = @(y) real(vec2im(A_rd'*vec(invF1(M(F1(y))))));

%% Running the iterations for each sigma
% Adding Guassian noise to the sinograms
y_rd_pad_noise = GaussianNoise( y_rd_pad, Sigma );
y_rd_noise     = GaussianNoise( y_rd,     Sigma );

%% Decimating y_rd_pad and theta_pad angles
y_rd_noise      = y_rd_noise(:,1:DecFactor:end);
y_rd_pad_noise  = y_rd_pad_noise(:,1:DecFactor:end);

[PSNR,SNR] = psnr(y_rd_noise,y_rd(:,1:DecFactor:end));

%% PP Interpolation
% Interpolating for each of the interpolation methods
for i=1:length(InterpMethod)
    y_pp{i} = InterpolateSinogram( y_rd_pad_noise, N, theta_pad, InterpMethod{i});
end
%%
tic;
y_fast = InterpolateFast(y_rd_pad_noise,N,theta_pad);
toc;