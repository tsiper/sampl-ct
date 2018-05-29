%% Resampling Radon to PP using our specific sub-space Kernel
% Comparing between different window sizes etc.
% Initializing:
Initialize;

%% PARAMETERS
N  = 28;                  % The image side resolution

% Chossing a phantom
PhantomTypes = {'brain','thorax','shepp','zubal'};
PhantomType  = PhantomTypes{1};

WindowSize =  [3,4,6,8,10,15];       % The radius of the window
% Sigma      = 0.05;    % Noise variance
SNR = 10:5:60;
Sigma = fliplr(10.^(-SNR/20));

Mp         = 2*N + 2;  % Number of sensors

% The default parameters for the kernel
% B         = floor(Mp/8);
B         = floor(Mp/6);
% R         = 1/pi/2/sqrt(2);
R         = 1/pi/5;  
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

%% PP Interpolation for every noise level
PSNR_Vals = zeros(length(Sigma),length(WindowSize));
SNR       = zeros(length(Sigma),1);

wb = MyWaitbar(0,'Calculating SNR Graph for Window Compare');
y_pp    = cell(length(Sigma),length(WindowSize));
y_pp_eq = y_pp;
for i=1:length(Sigma)
    % Adding noise
    y_rd_pad_n = GaussianNoise( y_rd_pad,Sigma(i) );
    y_rd_noise = GaussianNoise( y_rd,    Sigma(i) );

    % Effective SNR of the input sinogram calculation
    SNR(i) = psnr(trimcols(y_rd_noise,N),trimcols(y_rd,N));

    % Interpolating for all methods
    for j=1:length(WindowSize)
        y_pp{i,j} = InterpolateWindow( y_rd_pad_n, N, theta_pad, WindowSize(j), B, R, W );
        y_pp_eq{i,j} = EqualizeImage(y_pp{i,j},y_pp0);
        PSNR_Vals(i,j) = psnr(y_pp_eq{i,j},y_pp0);
    end
    
    % Updating the waitbar
    MyWaitbar(i/length(Sigma),wb); drawnow;    
end
close(wb);

%% Saving Results
save './FigureScripts/Param_Compare.mat' -regexp ^(?!(y_pp|y_pp_eq)$).;

%% loading the data from the saved file
load './FigureScripts/Param_Compare.mat'

%% Plotting the Param_Compare
WindowLegend = cell(1,length(WindowSize));
for i=1:length(WindowSize); WindowLegend{i}=['$K = ',num2str(WindowSize(i)),'$']; end

% %Plotting all
% figure;
% plot(SNR,PSNR_Vals);
% legend(WindowLegend);


SizeX = 1440; SizeY = 600;
TextParams = {'FontSize',28}; 
offset = 0;
figure('Position',[50 50 SizeX,SizeY],'Name','SNR_SNR_Compare');
LineSpec = {'+-.','o-','*--','^-','d-','s--'};
ColorSpec = {[.8 0 0],[0 0.5 0],[.7 0 .7],[0 0 0],[0 0 .8],[0 0.5 0.5]};
for i=1:length(LineSpec)
    h1 = subplot(121);plot(SNR,PSNR_Vals(:,i+offset),LineSpec{i},'LineWidth',2,'MarkerSize',10,'Color',ColorSpec{i});
    xlabel('Polar Sinogram PSNR'); ylabel('Resampled PP Sinogram PSNR');
    hold on;set(gca,TextParams{:});grid on;
    xlim([10 55]);
    rectangle('Position',[27,33,6,6], 'LineWidth',3);
    text(34,35,'Zoom Region','FontSize',26);

    h2 = subplot(122);ax2 = plot(SNR,PSNR_Vals(:,i),LineSpec{i},'LineWidth',2.5,'MarkerSize',10,'Color',ColorSpec{i});
    xlim([27 33]);ylim([32 38]);
    xlabel('Zoom-In PSNR'); ylabel('Zoom-In PP Sinogram PSNR');
    hold on;set(gca,TextParams{:});grid on;
end

l1 = legend(h1,WindowLegend{1+offset:length(LineSpec)+offset},'Location','southeast');
l2 = legend(h2,WindowLegend{1+offset:length(LineSpec)+offset},'Location','southeast');

set(l1,'Interpreter','latex')
set(l2,'Interpreter','latex')

%% Saving
SaveFigure('Window_Size_Compare_SNR',SizeX,SizeY);

