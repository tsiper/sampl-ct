%% Running a DEMO for SPURS resampling from a Polar trajectory

clear;
close all;
% clc;

setenv('SPURS_DIR', pwd);
setenv('SPURS_RUN_TIME', datestr(now,'ddmmyyyyTHHMMSS'));

%%
sqrtN = 35; % points per row / column of the cartesian grid

%==>phatom:
% PhantomString = 'Phantom1';
% PhantomString = 'AnalyticalSL';
PhantomString = 'Brain';

%==>Desired level of input noise in dB
% DesiredSNR = inf; % no noise
% DesiredSNR = 30;
Sigma = 0.015;
% DesiredSNR = 60;

%==>For a radial Trajectory:
Mproj = sqrtN*2;%The number of spokes of a radial trajectory, M=512*NSpokes.
Nray  = 2*ceil(norm([sqrtN sqrtN]-floor(([sqrtN sqrtN]-1)/2)-1))+3;
M = Mproj*Nray;

%==>SPURS Parameters
BsplineDegree = 3;
OverGridFactor = 2;
Rho = 1e-3;
Niterations = 5;
FilterInImageSpace = 1;

%% Genereate sampling and reconstruction grids

N = sqrtN^2;
% [ SamplingGridCoordinates ] = ConstructGridRadial0( Nbins, NSpokes);
PolarGrid = BuildPolarGrid(Nray,Mproj)/Nray*sqrtN;

%% Calculate phantom samples

Theta = (0:1/Mproj:1-1/Mproj)*180;

% [RefImg, PhantomSamples] = SamplePhantom(PhantomString,TrajrctoryString,sqrtN,M,SamplingGridCoordinates,0);
RefImg   = LoadPhantom(sqrtN,'zubal');
Sinogram = radon(RefImg,Theta);
Sinogram_Noisy = Sinogram + Sigma*randn(size(Sinogram))*range(Sinogram(:));
PhantomSamples = vec(F1(Sinogram_Noisy))/N;


%% Add noise to the samples

% OriginalPhantomSamples = PhantomSamples;
% [PhantomSamples,achivedSNR] = AddNoiseToData(PhantomSamples,DesiredSNR);

%% Run SPURS

SPURS_settings.sqrtN = sqrtN;
SPURS_settings.KernelFunctionString = 'Bspline';
SPURS_settings.KernelFunctionDegree = BsplineDegree;
SPURS_settings.ReusePrecalculatedData = 1;
SPURS_settings.Rho = Rho;
SPURS_settings.Niterations = Niterations;
SPURS_settings.UseW = 0;
SPURS_settings.ForceGenrateNewPhi = 0;
SPURS_settings.ForceFactorPsi = 0;
SPURS_settings.SavePSI = 0;
SPURS_settings.OverGridFactor = OverGridFactor;
SPURS_settings.alpha = 1;
SPURS_settings.CalcOptimalAlpha = 1;
SPURS_settings.FilterInImageSpace = FilterInImageSpace;

% tic
[ OutputImages, ReconstructedPhantomSamples] = SPURS(PhantomSamples, PolarGrid, SPURS_settings);
% toc
% disp(['SPURS using B-spline(',num2str(SPURS_settings.KernelFunctionDegree),') with Rho=',num2str(SPURS_settings.Rho),' finished. N=',num2str(N),' ,M=',num2str(M),', ISNR=',num2str(achivedSNR),' dB']);

PSNR_PreEqualizer = psnr(OutputImages(:,:,end),RefImg);

%% Equalizing the best image
ReconImg = zeros(sqrtN);
ReconImg(1:sqrtN-2,1:sqrtN-1) = OutputImages(3:end,2:end,end);
ReconImg(ReconImg<0) = 0;
ReconImg = ReconImg*fminsearch(@(A) -psnr(A*ReconImg,RefImg),1);

PSNR_Final = psnr(ReconImg,RefImg);

CompareImages(ReconImg,RefImg);

%% Analyze results
% TitleString = ['SPURS with N=',num2str(N),' ,M=',num2str(M),', ISNR=',num2str(achivedSNR),' dB'];
% SPURS_IQM = AnalyzeResult(OutputImages,RefImg);