ReconstructionParams.toOM=1;
%% Specs
ReconstructionParams.PhantomRes = 64;
save_name = 'test1';
ReconstructionParams.toRecord = 0;
ReconstructionParams.fileRecordName = 'new_moreCircles128_noOptim';
%% Spectra Initialization
SpectralParams.K =4;
SpectralParams.N0 = 10^20;
SpectralParams.Vtube = [120];
SpectralParams.Threshold = [30,40,65,100];
SpectralParams.InterpolationRes = 1000;
SpectralParams.ToOptimizeThresholds = 0;
OM_Vtube = 100;
%% choosing materials
% materials for MusFromMaterials func (will be inside BodyMaterials):
% adipose, blood, bone_compact, bone_cortical, brain, lung, muscle_skeletal,
% muscle_striated, skin, soft_tissue, water
I = 3;
ChosenMaterials = { 'I','Ca', 'soft_tissue'}; %choose the materials you wish for the basis
ReconstructionParams.base = [0.125 1 3]./255;%Weightening for the material densities


%% Loading phantom
% chosen phantom can be either: 'circles', 'more_circles', 'paint',
% 'spectral_shepp', 'spectral_zubal'.
chosen_phantom = 'more_circles';

%% Creating spectra
scan_type = 'Thresholding';

%% Spectral Reconstruction - using FISTA
global DebugFlag PlotFlag;
DebugFlag=1;
PlotFlag=1;
%Parameters for reconstruction 0.0001
ReconstructionParams.noiseLevel= 0.03;
ReconstructionParams.lambda = 0.001; %lambda for TV step
ReconstructionParams.L =0.005*SpectralParams.K; %Lipschitz coeeficient. optimal 90 with calibration region , 800 for distinct materials,0.5adipose&water
ReconstructionParams.La_coeffs =[0.01,0.01,20];% [1,0.001,500]; %a metric coefficients, used for scaling
ReconstructionParams.iters = 1000; %Maximal number of iterations of the minimization algorithm
ReconstructionParams.TV_Iters = 10; %Number of iterations in TV (broken code)
ReconstructionParams.L_coeffs = ones(1,SpectralParams.K); %P metric coefficients, used for scaling
ReconstructionParams.toTV =1; %Flag stating whether to perform TV or not

