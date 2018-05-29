%% Resampling Radon to PP using our specific sub-space Kernel

%% Initializing
Initialize;
DebugFlag = 0;
SaveFlag = 0;
tic;

%% PARAMETERS

% Chossing a phantom
PhantomTypes = {'brain','thorax','shepp','zubal'};
PhantomType  = PhantomTypes{4};

% The interpolation methods to run
% InterpMethod = {'filt'};
InterpMethod = {'filt'};
% InterpMethod = {'nearest','linear','pchip','spline','filt'};
% InterpMethod = {'linear','spline','filt'};

% Main run paramters
N          = 256;            % The image side resolution
WindowSize =  16;            % The radius of the window
% SNR = 30;              % Signal to noise ratio
% DesiredSNR = 25:3:65;
% DesiredSNR = 30:5:35;
DesiredSNR = 32;
Sigma      = 10.^(-DesiredSNR/20); % Noise variance
% Sigma      = 0;             % Noise variance
Mp         = 2*N + 2;       % Number of projection angles
DecFactor  = 3;             % Decimation factor of the measurements
% ScanType   = 'MATLAB';      % The scanning approach {'MATLAB' / 'AIR'}

% TV Fista Params
% PP Params
% lambda=6e-3; % Good for 128
lambda=8e-2;
% lambda = 0.05;
L = 6;
FISTA_Iters = 12;
TV_Iters = 25;
IterFactor = 8;


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
x_art_d = x_fbp; x_art_tv_d = x_fbp; x_rd_tv_d = x_fbp;
x_spurs = x_fbp; x_spurs_tv = x_fbp; x_fbp_tv = x_fbp;
y_pp   = cell(length(InterpMethod),length(Sigma));
x_pp   = y_pp; y_pp_d = y_pp; x_pp_d = y_pp; 
x_tv_d = y_pp; x_tv   = y_pp; 

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
At   = @(y) real(vec2im(A_rd'*vec(y)));

%% Running the iterations for each sigma
% wb = MyWaitbar(0,'Reconstructing all combinations');
for n=1:length(Sigma)
    % Adding Guassian noise to the sinograms
    y_rd_pad_noise = GaussianNoise( y_rd_pad, Sigma(n) );
    y_rd_noise     = GaussianNoise( y_rd,     Sigma(n) );
    

    %% Decimating y_rd_pad and theta_pad angles
    y_rd_noise      = y_rd_noise(:,1:DecFactor:end);
    
    y_rd_denoise   = DenoiseSinogram(y_rd_noise);

    y_rd_pad_noise  = y_rd_pad_noise(:,1:DecFactor:end); 
    
    [PSNR(n),SNR(n)] = psnr(y_rd_noise,y_rd(:,1:DecFactor:end));
    
    %% PP Interpolation
    % Interpolating for each of the interpolation methods
    for i=1:length(InterpMethod)
        y_pp{i,n} = InterpolateSinogram( y_rd_pad_noise, N, theta_pad, InterpMethod{i});
    end


    %% Computing Results
    Results(2*length(InterpMethod)+1,1,n) = psnr(y_rd_noise,y_rd(:,1:DecFactor:end));   % PSNR for FBP
    Results(2*length(InterpMethod)+1,2,n) = ssim(y_rd_noise/max(y_rd(:)),y_rd(:,1:DecFactor:end)/max(y_rd(:)));   % SSIM for FBP
    % Equalizing the results
    for  i=1:length(InterpMethod)
        [y_pp{i,n},ParamTemp] = EqualizeImage(y_pp{i,n},y_pp0);
        Results(i,1,n) = ParamTemp.PSNR;
        Results(i,2,n) = ParamTemp.SSIM;
    end

    %% Reconstructions
%     % Using Conjugate-Gradients to solve each of the Pseudo-Polar sinograms
%     for i=1:length(InterpMethod)
%         x_pp{i,n}   = Proj_C(real(conjgrad(@(x) App_T(M(App(x))), App_T(M(F1(y_pp{i,n}))) , zeros(N) ,cg_iters )));
% %         x_pp_d{i} = Proj_C(conjgrad(@(x) App_T(M(App(x))), App_T(M(F1(y_pp_d{i}))) , zeros(N) ,cg_iters ));
%     end

%     % Equalizing and saving quality metrics
%     for i=1:length(InterpMethod)
%         [x_pp{i,n}, TempParams] = EqualizeImage(x_pp{i,n},x0);
%         Results(i,3,n)     = TempParams.PSNR;
%         Results(i,4,n)     = TempParams.SSIM;
%     end

    x_fbp{n}    = Proj_C(vec2im(fbp(A_rd,y_rd_noise(:),theta)));
    x_fbp{n}    = imresize(x_fbp{n},(N+1)/N);
    x_fbp{n}    = x_fbp{n}(1:end-1,1:end-1);
    
    x_fbp_d{n}    = Proj_C(vec2im(fbp(A_rd,y_rd_denoise(:),theta)));
    x_fbp_d{n}    = imresize(x_fbp_d{n},(N+1)/N);
    x_fbp_d{n}    = x_fbp_d{n}(1:end-1,1:end-1);
    
%     end
    % Equalizing the image
    [x_fbp{n}, fbpParams] = EqualizeImage(ShiftImage(x_fbp{n},x0),x0);
    [x_fbp_d{n}, fbpParams] = EqualizeImage(ShiftImage(x_fbp_d{n},x0),x0);
    Results(2*length(InterpMethod)+1,3,n) = fbpParams.PSNR;
    Results(2*length(InterpMethod)+1,4,n) = fbpParams.SSIM;

    %% Starting SPURS
    % Fixing the angle range of the sinogram so it's from 0 to 180
%     y_fix_spurs = [y_rd_noise(:,N/2+2:end), flipud(y_rd_noise(:,1:N/2+1))];
    % y_fix = y_rd;
    PhantomSamples = vec(F1(y_rd_noise))/numel(y_rd_noise);
    % Fixing the bug in F1 (don't know why it happens
    % Loading default settings
    SPURS_settings = SPURS_DefaultParams(N);
    PolarGrid = BuildPolarGrid(N_rd,theta)/N_rd*N;
    [ OutputImages, ~] = SPURS(PhantomSamples, PolarGrid, SPURS_settings);
    
    
    PhantomSamples_d = vec(F1(y_rd_denoise))/numel(y_rd_denoise);
    % Loading default settings
    SPURS_settings = SPURS_DefaultParams(N);
    PolarGrid = BuildPolarGrid(N_rd,theta)/N_rd*N;
    [ OutputImages_d, ~] = SPURS(PhantomSamples_d, PolarGrid, SPURS_settings);
    
    
    
    %% SPURS Analysis
    [x_spurs{n},ParamSPURS] = EqualizeImage(ShiftImage(OutputImages(:,:,1),x0),x0);
    
    [x_spurs_d{n},ParamSPURS_d] = EqualizeImage(ShiftImage(OutputImages_d(:,:,1),x0),x0);
    
    Results(2*length(InterpMethod)+2,3,n) = ParamSPURS.PSNR;
    Results(2*length(InterpMethod)+2,4,n) = ParamSPURS.SSIM;


    %% Applying TV solvers to data
    D = DecOperator(DecFactor,'uniform');
%     A_Fista  = @(x) (D(Mnew(App(x,Mp/2))));
%     At_Fista = @(y) DecFactor*App_T(D(y));
    A_Fista  = @(x) D(M(App(x,Mp/2)));
    At_Fista = @(y) real(App_T(y));
    % x_tv  =  MFISTA(y_Fista,A_Fista,At_Fista,x_Fista,lambda,L,mfista_iters,TV_Iters,x0);
    for i=1:length(InterpMethod)
        y_Fista   = M(F1(y_pp{i,n}));
        x_tv{i,n}   = FISTA_TV(y_Fista,A_Fista,At_Fista,zeros(size(x0)),lambda,L,FISTA_Iters,TV_Iters,x0);
        x_tv{i,n}   = EqualizeImage(real(x_tv{i,n}),x0);
    end

    %% Calculate more results
    for i=1:length(InterpMethod)
        Results(i,5,n) = psnr(x_tv{i,n},x0);
        Results(i,6,n) = ssim(x_tv{i,n},x0);
%         Results(length(InterpMethod)+i,5,n) = psnr(x_tv_d{i},x0);
%         Results(length(InterpMethod)+i,6,n) = ssim(x_tv_d{i}/max(x0(:)),x0/max(x0(:)));
    end

    %% TV Denoising for SPURS and FBP
    % TV denoising for filtered back-projection
    x_fbp_tv{n} = EqualizeImage(FGP(x_fbp{n},lambda,TV_Iters),x0);
    Results(2*length(InterpMethod)+1,5,n) = psnr(x_fbp_tv{n},x0);
    Results(2*length(InterpMethod)+1,6,n) = ssim(x_fbp_tv{n},x0);

    x_fbp_tv_d{n} = EqualizeImage(FGP(x_fbp_d{n},lambda,TV_Iters),x0);
    
    
    % TV denoising for SPURS
    x_spurs_tv{n} = EqualizeImage(FGP(x_spurs{n},lambda,TV_Iters),x0);
    Results(2*length(InterpMethod)+2,5,n) = psnr(x_spurs_tv{n},x0);
    Results(2*length(InterpMethod)+2,6,n) = ssim(x_spurs_tv{n},x0);
    
    x_spurs_tv_d{n} = EqualizeImage(FGP(x_spurs_d{n},lambda,TV_Iters),x0);
    
    %% Comparing to other state of the art methods

    % ART Reconstruction
    x_art{n} = vec2im(sart(A_rd,y_rd_noise(:),FISTA_Iters));
    x_art{n} = ShiftImage(imresize(x_art{n},(N+1)/N),x0);
    x_art{n} = EqualizeImage(x_art{n}(1:end-1,1:end-1),x0);
    
    % ART Reconstruction for denoised sinogram
    x_art_d{n} = vec2im(sart(A_rd,y_rd_denoise(:),FISTA_Iters));
    x_art_d{n} = ShiftImage(imresize(x_art_d{n},(N+1)/N),x0);
    x_art_d{n} = EqualizeImage(x_art_d{n}(1:end-1,1:end-1),x0);
    
    
    % ART + TV
    x_art_tv{n} = zeros(N);
    for k=1:FISTA_Iters*IterFactor
        x_art_tv{n} = vec2im(sart(A_rd,y_rd_noise(:),1,x_art_tv{n}(:)));
        x_art_tv{n} = FGP(x_art_tv{n},lambda/FISTA_Iters/IterFactor,TV_Iters);
    end
%     x_art_tv{n} = EqualizeImage(x_art_tv{n},x0);
    x_art_tv{n} = ShiftImage(imresize(x_art_tv{n},(N+1)/N),x0);
    x_art_tv{n} = EqualizeImage(x_art_tv{n}(1:end-1,1:end-1),x0);
    
    % ART + TV + denoised
    x_art_tv_d{n} = zeros(N);
    for k=1:FISTA_Iters*IterFactor
        x_art_tv_d{n} = vec2im(sart(A_rd,y_rd_denoise(:),1,x_art_tv_d{n}(:)));
        x_art_tv_d{n} = FGP(x_art_tv_d{n},lambda/FISTA_Iters/IterFactor,TV_Iters);
    end
%     x_art_tv{n} = EqualizeImage(x_art_tv{n},x0);
    x_art_tv_d{n} = ShiftImage(imresize(x_art_tv_d{n},(N+1)/N),x0);
    x_art_tv_d{n} = EqualizeImage(x_art_tv_d{n}(1:end-1,1:end-1),x0);
    
    
    %% Running for the iterative Polar+Tv    
    L_radon = PowerMethod(@(x) At(A(x)),256);
 
    x_rd_tv{n} = FISTA_TV(y_rd_noise,A,At,zeros(size(x0)),lambda,L_radon,IterFactor*FISTA_Iters,TV_Iters,x0);
%     x_rd_tv{n} = EqualizeImage(x_rd_tv{n},x0);
    x_rd_tv{n} = ShiftImage(imresize(x_rd_tv{n},(N+1)/N),x0);
    x_rd_tv{n} = EqualizeImage(x_rd_tv{n}(1:end-1,1:end-1),x0);
    
    % Now for denoised sinogram
    x_rd_tv_d{n} = FISTA_TV(y_rd_denoise,A,At,zeros(size(x0)),lambda,L_radon,IterFactor*FISTA_Iters,TV_Iters,x0);
%     x_rd_tv{n} = EqualizeImage(x_rd_tv{n},x0);
    x_rd_tv_d{n} = ShiftImage(imresize(x_rd_tv_d{n},(N+1)/N),x0);
    x_rd_tv_d{n} = EqualizeImage(x_rd_tv_d{n}(1:end-1,1:end-1),x0);
    
    %% State-of-the-art Results
    Results(2*length(InterpMethod)+3,5,n) = psnr(x_art{n},x0);
    Results(2*length(InterpMethod)+3,6,n) = ssim(x_art{n},x0);
    Results(2*length(InterpMethod)+4,5,n) = psnr(x_art_tv{n},x0);
    Results(2*length(InterpMethod)+4,6,n) = ssim(x_art_tv{n},x0);
    Results(2*length(InterpMethod)+5,5,n) = psnr(x_rd_tv{n},x0);
    Results(2*length(InterpMethod)+5,6,n) = ssim(x_rd_tv{n},x0);
    
    %% Putting results in a nice table :)
    T{n} = array2table([(1:length(ResultsRows))',Results(:,:,n)],...
        'VariableNames',[{'ind'},ResultsCols],'RowNames',ResultsRows); %#ok<SAGROW>
    disp(sortrows(T{n},6));
    
    %%
%     MyWaitbar(n/length(Sigma),wb); drawnow;
end
% close(wb);
toc;
%% Saving
% save ./FigureScripts/Reconstruction_Compare.mat T Results SNR N DecFactor WindowSize

%% Loading

% load ./FigureScripts/Reconstruction_Compare.mat


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting Constants and Configurations
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PlotStyle = {'x-.','o-','^-','*--','d-','s--'};
PlotColor = {[.8 0 0],[0 0.5 0],[.7 0 .7],[0 0 0],[0 0 .8],[0 0.5 0.5]};
LineWidth = 2;
FontSize  = 28;
PlotParams = {'LineWidth',2.5,'MarkerSize',12}; 

% The values that we wanna plot
yPSNR_Vals = squeeze(Results(:,1,:))';
ySSIM_Vals = squeeze(Results(:,2,:))';
xPSNR_Vals = squeeze(Results(:,5,:))';
xSSIM_Vals = squeeze(Results(:,6,:))';
if length(SNR) > 10
    SNR_Range = 2:length(SNR)-2;
else
    SNR_Range  = 1:length(SNR);
end

ResultsRowsStrings = {'Nearest Neighbour','Linear','Cubic','Spline','Subspace',...
    'Denoised Nearest Neighbour','Denoised Linear','Denoised Cubic','Denoised Spline',...
    'Subspace','FBP+TV','SPURS+TV','ART','ART+TV','Polar+TV'};


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Displaying the different reconstructions with zoom for specific SNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 1;
Width = 600; Height = 900;
figure('Position',[50 50 Width Height]);
g = 0.01; h = 0; w = 0;

set(0,'defaulttextinterpreter','none');
ax(1) = subtightplot(3,2,1,g,h,w);imshowzoom(x0,'a');
ax(2) = subtightplot(3,2,2,g,h,w);imshowzoom(x_fbp_tv{n},'b');
ax(3) = subtightplot(3,2,3,g,h,w);imshowzoom(x_rd_tv{n},'c');
ax(4) = subtightplot(3,2,4,g,h,w);imshowzoom(x_art_tv{n},'d');
ax(5) = subtightplot(3,2,5,g,h,w);imshowzoom(x_spurs_tv{n},'e');
ax(6) = subtightplot(3,2,6,g,h,w);imshowzoom(x_tv{1,n},'f');
linkaxes(ax);
set(0,'defaulttextinterpreter','latex');

SaveFigure(['Reconstructions_D',num2str(DecFactor),'_',num2str(N),...
    '_',num2str(DesiredSNR(n)),'_F',num2str(IterFactor)],Width,Height);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Displaying the different reconstructions after denoising the measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 1;
Width = 600; Height = 900;
figure('Position',[50 50 Width Height]);
g = 0.01; h = 0; w = 0;

set(0,'defaulttextinterpreter','none');
ax(1) = subtightplot(3,2,1,g,h,w);imshowzoom(x0,'a');
ax(2) = subtightplot(3,2,2,g,h,w);imshowzoom(x_fbp_tv_d{n},'b');
ax(3) = subtightplot(3,2,3,g,h,w);imshowzoom(x_rd_tv_d{n},'c');
ax(4) = subtightplot(3,2,4,g,h,w);imshowzoom(x_art_tv_d{n},'d');
ax(5) = subtightplot(3,2,5,g,h,w);imshowzoom(x_spurs_tv_d{n},'e');
ax(6) = subtightplot(3,2,6,g,h,w);imshowzoom(x_tv{1,n},'f');
linkaxes(ax);
set(0,'defaulttextinterpreter','latex');

SaveFigure(['Reconstructions_Denoised_D',num2str(DecFactor),'_',num2str(N),...
    '_',num2str(DesiredSNR(n)),'_F',num2str(IterFactor)],Width,Height);

