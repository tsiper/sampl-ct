%% Specs
% The ground truth phantom side resolution (NxN phantom, N=PhantomRes)
ReconstructionParams.PhantomRes = 256;
% The name of the file to save the simulation result
save_name = 'test1';
ReconstructionParams.toRecord =0;
ReconstructionParams.fileRecordName = 'simRecord';
%% Spectra Initialization
% Number of spectra involved in each scan
SpectralParams.K = 1;
% Amount of photons emitted, proportional to dosage
SpectralParams.N0 = 10^2;
% The tube voltage, if we are on "switching" it should have length(K) elements
% else (thresholding) just the tube voltage of the single emitted spectrum
SpectralParams.Vtube = [70 140];
% The bin right limits markers, the first bin starts at 0
SpectralParams.Threshold = [64,80, 98 ,120];
% The analog energy interpolation factor for the integral approximation
SpectralParams.InterpolationRes = 1000;

SpectralParams.ToOptimizeThresholds = 0;

%% choosing materials
% materials for MusFromMaterials func (will be inside BodyMaterials):
% adipose, blood, bone_compact, bone_cortical, brain, lung, muscle_skeletal,
% muscle_striated, skin, soft_tissue, water

% materials for XrayMu func will be string with Chemical formula or simple
% atomic notation

I = 3; %Number of materials
ChosenMaterials = { 'I','Ca','soft_tissue'}; %choose the materials you wish for the basis
ReconstructionParams.base = [1 1 1]./255; %Weightening for the material densities

%% Loading phantom
% chosen phantom can be either: 'circles', 'more_circles', 'paint',
% 'spectral_shepp', 'spectral_zubal'.
chosen_phantom = 'more_circles';

%% Creating spectra
%Spectra:
%scantype defines the type of the scan for the simulation. It can be:
%'switching' - used mainly for DECT, where the tube's voltage changes (usually 2 different
% voltages). This option simulates the needed spectra.
%'Thresholding' - used for spectral scan, where the detectors are energy
%sensitive. This option simulates the thresholded bins' spectra.
%'ExtSpectra' - used when there is an external measured spectra, and want to
% simulate the corresponding attenuation (the spectrum is interpulated
% through GetSpectra()).
scan_type = 'switching';

%% Spectral Reconstruction - using FISTA
global DebugFlag PlotFlag;
DebugFlag=1;
PlotFlag=1;

%Parameters for reconstruction
ReconstructionParams.noiseLevel = 0.001;
ReconstructionParams.toStop = 1;
ReconstructionParams.lambda = 0.001; %lambda for TV step
ReconstructionParams.L =200*SpectralParams.K; %Lipschitz coeeficient. optimal 90 with calibration region , 800 for distinct materials,0.5adipose&water
ReconstructionParams.La_coeffs = [0.125 1];%[0.125 1 20]; %a metric coefficients, used for scaling
ReconstructionParams.iters = 1000; %Maximal number of iterations of the minimization algorithm
ReconstructionParams.TV_Iters = 20; %Number of iterations in TV (broken code)
ReconstructionParams.L_coeffs = ones(1,SpectralParams.K); %P metric coefficients, used for scaling
ReconstructionParams.toTV =1; %Flag stating whether to perform TV or not
