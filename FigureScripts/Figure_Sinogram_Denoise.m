%% Plots a comparison of the different sinogram denoising results

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
InterpMethod = {'nearest','linear','pchip','spline','filt'};
% InterpMethod = {'linear','spline','filt'};

% Main run paramters
N          = 128;            % The image side resolution
WindowSize =  8;            % The radius of the window
% SNR = 30;              % Signal to noise ratio
DesiredSNR = 25:3:55;
% DesiredSNR = 30:5:35;
% DesiredSNR = 35;
Sigma      = 10.^(-DesiredSNR/20); % Noise variance
% Sigma      = 0;             % Noise variance
Mp         = 2*N + 2;       % Number of projection angles
cg_iters   = 3;             % Maximum iterations for conjugate gradient alg.
DecFactor  = 1;             % Decimation factor of the measurements
% ScanType   = 'MATLAB';      % The scanning approach {'MATLAB' / 'AIR'}

% The default parameters for the kernel
% B         = (Mp./(3:3:30));
% B = 4*pi*[1:0.5:2,3,10,20];
B = 4*pi*[(1:0.25:10),11:30];
% R         = 1/pi/2/sqrt(2);
R         = (1/pi/2);
% R         = 1/pi/2;
W         = (pi)*(2*N+1);


%% Initializing Variables and cells
% Loading the Phantom
x0 = LoadPhantom(N,PhantomType);

% Calculating the ideal Pseudo-Polar sinogram
y_pp0 = real(invF1(App(x0)));

% The real SNR and PSNR vectors
SNR_out = zeros(length(Sigma),length(B));
SNR_in  = zeros(size(Sigma));

%% Generating the scanned data
% The angle ranges
N_rd      = 2*ceil(sqrt(2)*(N-floor((N-1)/2)-1))+3; % Number of sensors
theta     = (0:1/Mp:1-1/Mp)*180;             % Radon Angles
theta_pad = [theta,theta+180];


% The projections
y_rd_pad = radon(x0,theta_pad);
y_rd_pad = padcols(y_rd_pad,2*N+1);

y_rd     = y_rd_pad(:,Mp/2+1:3*Mp/2);


figure;
subplot(121);imagesc(y_rd_pad);
subplot(122);imagesc(y_rd);


%% Running the iterations for each sigma

for i=1:length(Sigma)
    % Adding Guassian noise to the sinograms
    y_rd_pad_noise{i} = GaussianNoise( y_rd_pad, Sigma(i) ); %#ok<*SAGROW>
    y_rd_noise{i}     = GaussianNoise( y_rd,     Sigma(i) );

    for j=1:length(B)
        [~,SNR_in(i)] = psnr(y_rd_noise{i},y_rd);
        y_d_pad{i,j}  = DenoiseSinogram(y_rd_pad_noise{i},B(j),R,W,WindowSize);
        y_d{i,j}      = y_d_pad{i,j}(:,Mp/2+1:3*Mp/2);
        SNR_out(i,j)  = psnr(y_d{i,j},y_rd);
    end
end

%% Plotting
Width = 1440;
Height  = 720;
gap = 0.05; mh = [0.08 0.02]; mv=[0.15 0.05];
FontSize = 26;
LineWidth = 2.5;
PlotStyle  = {'x-.','o-','^-','.--','s-','d--'};
PlotStyle2 = {'--','-.',':','--','.--','s--'};
PlotColor = {[.8 0 0],[0 0.5 0],[.7 0 .7],[0 0 0],[0 0 .8],[0 0.5 0.5]};


figure('Position',[50 50 Width Height]);
B_range = [3,5,9,21,37];
subtightplot(1,2,1,gap,mv,mh);
for k=1:length(B_range)
    plot(SNR_in,SNR_out(:,B_range(k)),PlotStyle{k},'LineWidth',LineWidth, ...
        'MarkerSize',10,'Color',PlotColor{k});
    hold on;
end
plot(SNR_in,SNR_in,'k','LineWidth',LineWidth);
legend(cellfun(@(x)['$B=$',num2str(x)],num2cell(B(B_range)/(8*pi)),'UniformOutput',0),'location','southeast');
xlim(round([min(SNR_in),max(SNR_in)])); ylim(round([min(SNR_in)+5,max(SNR_in)-5]));
xlabel('Sinogram SNR [dB]');ylabel('Denoised Sinogram SNR [dB]');
set(gca,'FontSize',FontSize);grid on;

SNR_Range = 2:2:10;
subtightplot(1,2,2,gap,mv,mh);
k = 1;
for SNR_Point=SNR_Range
    plot(B/8/pi,SNR_out(SNR_Point,:),PlotStyle2{k},'LineWidth',LineWidth,...
        'MarkerSize',10,'Color',PlotColor{k},'MarkerSize',12);
    k = k+1;
    hold on;
end
xlabel('Value of $B$'); grid on;
ylim(round([min(SNR_in)+5,max(SNR_in)-5]));xlim([0 10]);
legend(cellfun(@(x)['SNR ',num2str(x),'dB'],num2cell(round(SNR_in(SNR_Range))),'UniformOutput',0),'location','southeast');
set(gca,'FontSize',FontSize);

SaveFigure('B_Analysis',Width,Height);
