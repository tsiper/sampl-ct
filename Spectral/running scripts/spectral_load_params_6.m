%% Specs
ReconstructionParams.PhantomRes = 128;
ReconstructionParams.toRecord = 1;
ReconstructionParams.fileRecordName = 'run_6';
save_name = 'run_6';
%% Spectra Initialization
SpectralParams.K = 4;
SpectralParams.N0 = 10^20;
SpectralParams.Vtube = [120];
SpectralParams.Threshold = [64,80, 98 ,120];
SpectralParams.InterpolationRes = 1000;
SpectralParams.ToOptimizeThresholds = 1;

%% choosing materials
% materials for MusFromMaterials func (will be inside BodyMaterials):
% adipose, blood, bone_compact, bone_cortical, brain, lung, muscle_skeletal,
% muscle_striated, skin, soft_tissue, water

IRecon = 2; %Number of materials to reconstruct
IReal = 3; %Number of materials actually exist
ChosenMaterialsRecon = {'I', 'Ca'}; %choose the materials you wish to reconstruct
ChosenMaterialsReal = {'soft_tissue','I','Ca'}; %choose the materials of the phantom
ReconstructionParams.base = [0.125 1]./255;%Weightening for the reconstruction material densities
ReconstructionParams.baseReal = [3 0.125 1]./255;%Weightening for the real material densities

chosen_phantom = 'circles';
scan_type = 'Thresholding';

%% Spectral Reconstruction - using FISTA
global DebugFlag PlotFlag;
DebugFlag=1;
PlotFlag=1;

%Parameters for reconstruction

ReconstructionParams.noiseLevel = 0;%0.01;
ReconstructionParams.lambda = 0.000001; %lambda for TV step
ReconstructionParams.L =100*SpectralParams.K; %Lipschitz coeeficient. optimal 90 with calibration region , 800 for distinct materials,0.5adipose&water
ReconstructionParams.La_coeffs = [0.125 1]; %a metric coefficients, used for scaling
ReconstructionParams.iters = 1000; %Maximal number of iterations of the minimization algorithm
ReconstructionParams.TV_Iters = 20; %Number of iterations in TV (broken code)
ReconstructionParams.L_coeffs = ones(1,SpectralParams.K); %P metric coefficients, used for scaling
ReconstructionParams.toTV =1; %Flag stating whether to perform TV or not
ReconstructionParams.toStop =0; %Flag stating whether to stop with stopping criterion or not
