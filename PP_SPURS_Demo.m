 % clear all;
% close all;
Initialize;
rmpath(genpath([pwd,'/Analysis/']));
% clc;


%% SHAHAR - MY VARIABLES

% SamplingGridCoordinates [x,y] coordinate pairs
% 
%%

[workdirpath, name, ext] = fileparts(mfilename('fullpath'));
workdirpath = [workdirpath,'\Toolboxes\Spurs'];
% cd(workdirpath);
setenv('SPURS_DIR', workdirpath);
setenv('SPURS_RUN_TIME', datestr(now,'ddmmyyyyTHHMMSS'));

addpath(genpath(workdirpath));

%% User defined settings and parameters: 

sqrtN = 128; % points per row / column of the cartesian grid
DecFactor = 6;
N = sqrtN^2;

%==>phatom:
% PhantomString = 'Phantom1';
% PhantomString = 'AnalyticalSL';
PhantomString = 'Brain';

%==>Desired level of input noise in dB
% DesiredSNR = inf; % no noise
DesiredSNR = 120;
% DesiredSNR = 60;

%==>Trajectory:
% TrajrctoryString = 'Radial';
TrajectoryString = 'Spiral';

%==>For a radial Trajectory:
% NSpokes = 60;%The number of spokes of a radial trajectory, M=512*NSpokes.
% NSpokes = 120;
Nbins = sqrtN*2;

%==>For a spiral Trajectory:
% M = 20000; %The nember of sampling points for a spiral trajectory.
% M = 30000;
% M = 50000;

%==>SPURS Parameters
BsplineDegree = 4;
OverGridFactor = 2;
Rho = 1e-5;
Niterations = 5;
FilterInImageSpace = 1;

%% Building the PP grid and building its phantom

[PPGrid,theta_vec] = BuildPPGrid(sqrtN);

x = LoadPhantom(sqrtN,'zubal');
y = -App(flipud(x))/((sqrtN+1/2)*(sqrtN+1));

% Limited Angle Problem
% y = y(:,16:end);
% PPGrid = PPGrid(15*(2*sqrtN+1)+1:end,:);

% Weird scaling factor BUG TODO:Investigate
PPGrid = PPGrid/2;

[SamplingGridCoordinates,y_dec] = DecimateSPURS(PPGrid,y,DecFactor);

RefImg = x;
x0 = RefImg;
PhantomSamples = [SamplingGridCoordinates(:,1), SamplingGridCoordinates(:,2), real(y_dec(:)), imag(y_dec(:))];

M = size(SamplingGridCoordinates,1);
Ns = sqrt(M/pi());

%% Genereate sampling and reconstruction grids

% N = sqrtN^2;
% switch (TrajrctoryString)
%     case 'Radial'
%         [ SamplingGridCoordinates ] = ConstructGridRadial0( Nbins, NSpokes);
%     case 'Spiral'
%         [ SamplingGridCoordinates ]  = ConstructGrid1ArmSpiral( sqrtN, M, Ns);
% end;

% Correct OverGridFactor value
if OverGridFactor ~= ceil(sqrtN*OverGridFactor/2)*2/sqrtN
    disp(['Over grid factor was corrected from ',num2str(OverGridFactor),' to ',num2str(ceil(sqrtN*OverGridFactor/2)*2/sqrtN)]);
    OverGridFactor = ceil(sqrtN*OverGridFactor/2)*2/sqrtN;
end

% ReconstructionGridCoordinates = ConstructGridCartesian( sqrtN, 1, 0, 0, OverGridFactor);
ReconstructionGridCoordinates = 2*ConstructGridCartesian( sqrtN, 1, 0, 0, OverGridFactor);


%% Calculate phantom samples

% [RefImg, PhantomSamples] = SamplePhantom(PhantomString,TrajrctoryString,sqrtN,M,SamplingGridCoordinates,0);

% ImageGrid = (-sqrtN/2:sqrtN/2-1);
ImageGrid = (-sqrtN:sqrtN);
% figure(2);

%% Add noise to the samples

OriginalPhantomSamples = PhantomSamples;
[PhantomSamples,achivedSNR] = AddNoiseToData(PhantomSamples,DesiredSNR);

%% Run SPURS

SPURS_settings.sqrtN = sqrtN;
SPURS_settings.KernelFunctionString = 'Bspline';
SPURS_settings.KernelFunctionDegree = BsplineDegree;
SPURS_settings.ReusePrecalculatedData = 0;
SPURS_settings.Rho = Rho;
SPURS_settings.Niterations = Niterations;
SPURS_settings.UseW = 0;
SPURS_settings.ForceGenrateNewPhi = 1;
SPURS_settings.ForceFactorPsi = 0;
SPURS_settings.SavePSI = 0;
SPURS_settings.OverGridFactor = OverGridFactor;
SPURS_settings.alpha = 1;
SPURS_settings.CalcOptimalAlpha = 1;
SPURS_settings.FilterInImageSpace = FilterInImageSpace;

b = PhantomSamples(:,3)+1i*PhantomSamples(:,4);
% tic
[ OutputImages, ReconstructedPhantomSamples] = SPURS(b, SamplingGridCoordinates, SPURS_settings);
% toc
disp(['SPURS using B-spline(',num2str(SPURS_settings.KernelFunctionDegree),') with Rho=',num2str(SPURS_settings.Rho),' finished. N=',num2str(N),' ,M=',num2str(M),', ISNR=',num2str(achivedSNR),' dB']);

%% Analyze results
TitleString = ['SPURS with N=',num2str(N),' ,M=',num2str(M),', ISNR=',num2str(achivedSNR),' dB'];
SPURS_IQM = AnalyzeResult(OutputImages,RefImg,TitleString);


%% Plotting Stuff

% Plotting the grids
figure('Name','Sampling and reconstruction Grids','NumberTitle','off');
plot(SamplingGridCoordinates(:,1),SamplingGridCoordinates(:,2),'o');
title({['Sampling on a ',TrajectoryString,' trajectory with M= ',num2str(M)],['Reconstrucion on a cartesian grid with \sigma=',num2str(OverGridFactor),', N=',num2str(N)]});
grid on; grid minor
axis([-1/16 0 -1/16 0].*sqrtN); axis square;
hold on;
plot(ReconstructionGridCoordinates(:,1),ReconstructionGridCoordinates(:,2),'x');
legend('Sampling Point','Reconstruction Point');

% Showing the Phantom
figure('Name','Expected phantom Image','NumberTitle','off');
imagesc(ImageGrid,ImageGrid,RefImg);axis square; colormap('gray');
title({['Expected result'],[PhantomString,' Phantom (',num2str(sqrtN),'X',num2str(sqrtN),') image']});axis off;


%% PSNR Compare after cropping
x_hat = OutputImages(:,:,1);
x_hat(x_hat<0) = 0;
x_hat(x_hat>1) = 1;
figure;imagesc(x_hat);
psnr(x_hat,x0);