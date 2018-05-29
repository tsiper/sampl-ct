%% Resampling Radon to PP using our specific sub-space Kernel

%% Initializing
Initialize;

%% PARAMETERS
N  = 64;                  % The image side resolution

% Chossing a phantom
PhantomTypes = {'brain','thorax','shepp','zubal'};
PhantomType  = PhantomTypes{1};

WindowSize =  8;       % The radius of the window
% Sigma      = 0.05;    % Noise variance
Sigma = 0.025;

Mp         = 2*N + 2;  % Number of sensors
B         = floor(Mp/20);
R        = 1/pi/6;
W         = (pi)*(2*N+1);

%% Generating the scanned data
% Loading the Phantom
x0 = LoadPhantom(N,PhantomType);

% Number of sensors
N_rd     = 2*ceil(sqrt(2)*(N-floor((N-1)/2)-1))+3;

% Designing the scanning angles
theta      = (-1/4:1/Mp:3/4-1/Mp)*180;
theta_pad  = (-1/2:1/Mp:1-1/Mp)*180;

% Scanning using Matlab's radon function
y_rd     = radon(x0,theta);
y_rd_pad = padcols( radon(x0,theta_pad), 2*N+1 );

% The "true" Pseudo-Polar sinogram, calculated using the operator 'App'
y_pp0 = real(invF1(App(x0)));

%%Adding noise
y_rd_pad_n = GaussianNoise(y_rd_pad,Sigma);

%% Now for the Denoising

y_denoise = DenoiseSinogram(y_rd_pad_n);

%%
[ ppSinogram1 ,PSNR_Vals1] = InterpolateAll( y_denoise, N, theta_pad, WindowSize, B, R, W, y_pp0 );
[ ppSinogram2 ,PSNR_Vals2] = InterpolateAll( y_rd_pad_n, N, theta_pad, WindowSize, B, R, W, y_pp0 );
InterpMethods = {'nearest','linear','pchip','spline','subspace'};

%% Plotting
CompareImages(y_denoise,y_rd_pad);
CompareImages(y_rd_pad_n,y_rd_pad);

disp([PSNR_Vals1,PSNR_Vals2]);