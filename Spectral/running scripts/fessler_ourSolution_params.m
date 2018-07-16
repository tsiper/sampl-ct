%% Specs
% The ground truth phantom side resolution (NxN phantom, N=PhantomRes)
ReconstructionParams.PhantomRes = 256;
% The name of the file to save the simulation result
ReconstructionParams.toRecord = 1;
ReconstructionParams.fileRecordName = 'long_fessler_v5';
save_name = 'long_fessler_v5';
%% Spectra Initialization
% Number of spectra involved in each scan
SpectralParams.K = 4;
% Amount of photons emitted, proportional to dosage
SpectralParams.N0 = 2e20;
% The tube voltage, if we are on "switching" it should have length(K) elements
% else (thresholding) just the tube voltage of the single emitted spectrum
SpectralParams.Vtube = 140;
% The bin right limits markers, the first bin starts at 0
SpectralParams.Threshold = [50,51,52,53,54 , 70, 110,140];
% The analog energy interpolation factor for the integral approximation
SpectralParams.InterpolationRes = 1000;

SpectralParams.ToOptimizeThresholds = 1;

%% choosing materials
% materials for MusFromMaterials func (will be inside BodyMaterials):
% adipose, blood, bone_compact, bone_cortical, brain, lung, muscle_skeletal,
% muscle_striated, skin, soft_tissue, water

% materials for XrayMu func will be string with Chemical formula or simple
% atomic notation

% long&fessler: fat, blood, omnipaque300 (io-dine-based contrast agent), cortical bone,
% and air.
% 70keV "mono"
% resolution: 512x512
% GE LightSpeed X-ray CT fan-beam system
% switching: 140,80 keV, (normalized to get same total dose)
% # of photons: 2e5 in 140kV, 6e4 in 80kV

% For this simulation we let the triplet material library con-tainfive triplets
%selected fromfi ve materials: fat, blood, omni-paque300, cortical bone, and air,
%excluding the combination of omnipaque300 and cortical bone and the combination
%of omni-paque300 and fat.

%omnipaque 300: 300mg per 1mL => 0.3 g/cm3 iodine, and the rest is water.


I = 4; %Number of materials
ChosenMaterials = { 'adipose','blood','I','bone_cortical'};%,'air'}; %choose the materials you wish for the basis
% base is densities:
ReconstructionParams.base = [0.9,1.06,0.3,1.92]/ReconstructionParams.PhantomRes*64;%,0];
%Weightening for the material densities


%% Loading phantom
% chosen phantom can be either: 'circles', 'more_circles', 'paint',
% 'spectral_shepp', 'spectral_zubal'.
chosen_phantom = 'spectral_NCAT_512';

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
scan_type = 'Thresholding';

%% Spectral Reconstruction - using FISTA
global DebugFlag PlotFlag;
DebugFlag=1;
PlotFlag=1;

%Parameters for reconstruction

ReconstructionParams.noiseLevel = 0;
ReconstructionParams.lambda = 1e-7; %lambda for TV step
ReconstructionParams.L =500*SpectralParams.K; %100 worked for 64px %Lipschitz coeeficient. optimal 90 with calibration region , 800 for distinct materials,0.5adipose&water
ReconstructionParams.La_coeffs = [8,8,0.8,0.6];%[8,20,0.8,1.2]; worked for 64 %a metric coefficients, used for scaling
ReconstructionParams.iters = 20000; %Maximal number of iterations of the minimization algorithm
ReconstructionParams.TV_Iters = 20; %Number of iterations in TV (broken code)
ReconstructionParams.L_coeffs = ones(1,SpectralParams.K); %P metric coefficients, used for scaling
ReconstructionParams.toTV =1; %Flag stating whether to perform TV or not
ReconstructionParams.toStop =0; %Flag stating whether to stop with stopping criterion or not