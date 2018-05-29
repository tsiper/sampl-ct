%% Spectral CT
% Project by Gil Ben Ari and Rotem Zamir
clc;
clear all;
Initialize;
% to run script you need to add all the CT folder to the path.

%% Loading simulation parameters 
% Specs
% The ground truth phantom side resolution (NxN phantom, N=PhantomRes)
ReconstructionParams.PhantomRes = 128;
% The name of the file to save the simulation result
ReconstructionParams.toRecord = 0;
ReconstructionParams.fileRecordName = 'Circles128_TV0-00005_Noised';
save_name = 'Circles128_TV0-00005_Noised';

% Spectra Initialization
% Number of spectra involved in each scan
SpectralParams.K = 4;
% Amount of photons emitted, proportional to dosage
SpectralParams.N0 = 10^20;
% The tube voltage, if we are on "switching" it should have length(K) elements
% else (thresholding) just the tube voltage of the single emitted spectrum
SpectralParams.Vtube = [120];
% The bin right limits markers, the first bin starts at 0
SpectralParams.Threshold = [64,80, 98 ,120];
% The analog energy interpolation factor for the integral approximation
SpectralParams.InterpolationRes = 1000;

SpectralParams.ToOptimizeThresholds = 1;

% choosing materials
% materials for MusFromMaterials func (will be inside BodyMaterials):
% adipose, blood, bone_compact, bone_cortical, brain, lung, muscle_skeletal,
% muscle_striated, skin, soft_tissue, water

% materials for XrayMu func will be string with Chemical formula or simple
% atomic notation

I = 3; %Number of materials
ChosenMaterials = { 'I','Ca', 'soft_tissue'}; %choose the materials you wish for the basis
ReconstructionParams.base = [0.125 1 3]./255;%Weightening for the material densities


% Loading phantom
% chosen phantom can be either: 'circles', 'more_circles', 'paint',
% 'spectral_shepp', 'spectral_zubal'.
chosen_phantom = 'circles';

% Creating spectra
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

% Spectral Reconstruction - using FISTA
global DebugFlag PlotFlag;
DebugFlag=1;
PlotFlag=1;

%Parameters for reconstruction

ReconstructionParams.noiseLevel = 0.05;
ReconstructionParams.lambda = 0.00001; %lambda for TV step
ReconstructionParams.L =100*SpectralParams.K; %Lipschitz coeeficient. optimal 90 with calibration region , 800 for distinct materials,0.5adipose&water
ReconstructionParams.La_coeffs = [0.125 1 20]; %a metric coefficients, used for scaling
ReconstructionParams.iters = 200; %Maximal number of iterations of the minimization algorithm
ReconstructionParams.TV_Iters = 20; %Number of iterations in TV (broken code)
ReconstructionParams.L_coeffs = ones(1,SpectralParams.K); %P metric coefficients, used for scaling
ReconstructionParams.toTV =1; %Flag stating whether to perform TV or not
ReconstructionParams.toStop =0; %Flag stating whether to stop with stopping criterion or not;

%% Loading phantom
RGBphantom = LoadPhantom(ReconstructionParams.PhantomRes,chosen_phantom);
phantom = cell(I,1);
for ii=1:I
    phantom{ii} = ReconstructionParams.base(ii)*RGBphantom(:,:,ii);
end

%% Creating spectra
% calling GetSpectra or loading external spectra data
Spectra = GetSpectra(scan_type, SpectralParams, ChosenMaterials, ReconstructionParams.base);

%% Creating the scan of the object (simulating a spectral scan)
[ Spectra, projections_clean, RadonMaterials, ChosenMaterials ] = scan( Spectra, phantom, ChosenMaterials );

%% Adding Noise
projections_Gaussian_1 = cell(SpectralParams.K,1);
projections_Gaussian_3 = cell(SpectralParams.K,1);
projections_Gaussian_5 = cell(SpectralParams.K,1);

for n=1:SpectralParams.K
	projections_Gaussian_1{n} = GaussianNoise(projections_clean{n},0.01);
	projections_Gaussian_3{n} = GaussianNoise(projections_clean{n},0.03);
	projections_Gaussian_5{n} = GaussianNoise(projections_clean{n},0.05);
end

%% Defining the Tensors from the scan cell data

% %Define Clean Projection Tensor
% Projection_Clean_Tensor = zeros([size(projections_clean{1}),length(projections_clean)]);
% for i=1:length(projections_clean)
%     Projection_Clean_Tensor(:,:,i) = projections_clean{i};
% end
% % Defining the Projection Tensor With the Tensor-ToolBox
% Tensor_Clean_Proj=tensor(Projection_Clean_Tensor); %% define the tensor with the toolbox
% 
% 
% %Define noised Projections Tensor - Gaussian Noise 0.01
% Projection_Gaussian_1_Tensor = zeros([size(projections_clean{1}),length(projections_clean)]);
% for i=1:length(projections_clean)
%     Projection_Gaussian_1_Tensor(:,:,i) = projections_Gaussian_1{i};
% end
% % Defining the Projection Tensor With the Tensor-ToolBox
% Tensor_Gaussian_1_Proj=tensor(Projection_Gaussian_1_Tensor); %% define the tensor with the toolbox
% 
% 
% %Define noised Projections Tensor - Gaussian Noise 0.03
% Projection_Gaussian_3_Tensor = zeros([size(projections_clean{1}),length(projections_clean)]);
% for i=1:length(projections_clean)
%     Projection_Gaussian_3_Tensor(:,:,i) = projections_Gaussian_3{i};
% end
% % Defining the Projection Tensor With the Tensor-ToolBox
% Tensor_Gaussian_3_Proj=tensor(Projection_Gaussian_3_Tensor); %% define the tensor with the toolbox
% 

%Define noised Projections Tensor - Gaussian Noise 0.05
Projection_Gaussian_5_Tensor = zeros([size(projections_clean{1}),length(projections_clean)]);
for i=1:length(projections_clean)
    Projection_Gaussian_5_Tensor(:,:,i) = projections_Gaussian_5{i};
end
% Defining the Projection Tensor With the Tensor-ToolBox
Tensor_Gaussian_5_Proj=tensor(Projection_Gaussian_5_Tensor); %% define the tensor with the toolbox

%%
number_of_ranks = 300;

%% Decomposition of the clean projection
% collect the data with many CP rank

TensorProj_Clean_CP_ALL = cell(1,number_of_ranks);
TensorProj_Clean_ALL = repmat({double(Tensor_Clean_Proj)},1,number_of_ranks);

% Starting the waitbar
counter = ConsoleProgressBar();
% Starting main iteration of MFISTA
counter.start(); counter.setText('Finding the best CP Decomposition Rank..');

for rank = 1:number_of_ranks
    TensorProj_Clean_CP_ALL{1,rank} = double(cp_als(Tensor_Clean_Proj,rank));
    % Collect data of the PSNR between the Tensor+Noise+CP and the clean Tensor
    PSNR_vec_Clean_ALL(rank) = psnr(TensorProj_Clean_CP_ALL{1,rank},TensorProj_Clean_ALL{1,rank});
    counter.setValue(rank/number_of_ranks);
end
counter.stop();

%find the max PSNR (value and index)
[maxPSNR_Clean_val_ALL , CP_Rank_Clean] = max(PSNR_vec_Clean_ALL(:));

%% Graph Section
%collect data for the comparing
%Calc the /PSNR/ of ALL proj
PSNR_ALL_Clean_NoCP = psnr(double(Tensor_Clean_Proj),double(Tensor_Clean_Proj));
% create /PSNR/ VECTOR (without CP) for the graph of ALL
PSNR_vec__Clean_NoCp_ALL = repmat(PSNR_ALL_Clean_NoCP,1,number_of_ranks);

% Create Graph of the Compare Between The Ground Truth and the Projections with several CP RANK.
figure(15);

%ALL projection
plot(1:number_of_ranks,PSNR_vec_Clean_ALL);
title('PSNR of the Proj after different decompositions');
xlim([0 number_of_ranks]);
xlabel('Rank');
ylabel('PSNR');
grid on;
hold on;

% mark the best RMSE/PSNR/SSIM  
plot(CP_Rank_Clean , maxPSNR_Clean_val_ALL,'--gs'); 
hold on;

%adding plots of Reference(ALL) (without CP)
plot(1:number_of_ranks,PSNR_vec__Clean_NoCp_ALL);
hold on;

legend('PSNR',['PSNR BEST= (',num2str(CP_Rank_Clean),' , ',num2str(maxPSNR_Clean_val_ALL),')'],...
   ['PSNR NoCp = ',num2str(PSNR_ALL_Clean_NoCP)]);
%% Decomposition of the noised (gaussian 0.01) projection
% collect the data with many CP rank

TensorProj_Gaussian_1_CP_ALL = cell(1,number_of_ranks);
TensorProj_Clean_ALL = repmat({double(Tensor_Clean_Proj)},1,number_of_ranks);

% Starting the waitbar
counter = ConsoleProgressBar();
% Starting main iteration of MFISTA
counter.start(); counter.setText('Finding the best CP Decomposition Rank..');

for rank = 1:number_of_ranks
    TensorProj_Gaussian_1_CP_ALL{1,rank} = double(cp_als(Tensor_Gaussian_1_Proj,rank));
    % Collect data of the PSNR between the Tensor+Noise+CP and the clean Tensor
    PSNR_vec_Gaussian_1_ALL(rank) = psnr(TensorProj_Gaussian_1_CP_ALL{1,rank},TensorProj_Clean_ALL{1,rank});
    counter.setValue(rank/number_of_ranks);
end
counter.stop();

%find the max PSNR (value and index)
[maxPSNR_Gaussian_1_val_ALL , CP_Rank_Gaussian_1] = max(PSNR_vec_Gaussian_1_ALL(:));


%% Graph Section
%collect data for the comparing
%Calc the /PSNR/ of ALL proj
PSNR_ALL_NoCP = psnr(double(Tensor_Gaussian_1_Proj),double(Tensor_Clean_Proj));
% create /PSNR/ VECTOR (without CP) for the graph of ALL
PSNR_vec_NoCp_ALL = repmat(PSNR_ALL_NoCP,1,number_of_ranks);

% Create Graph of the Compare Between The Ground Truth and the Projections with several CP RANK.
figure(16);

%ALL projection
plot(1:number_of_ranks,PSNR_vec_Gaussian_1_ALL);
title('PSNR of the Proj after different decompositions');
xlim([0 number_of_ranks]);
xlabel('Rank');
ylabel('PSNR');
grid on;
hold on;

% mark the best RMSE/PSNR/SSIM  
plot(CP_Rank_Gaussian_1 , maxPSNR_Gaussian_1_val_ALL,'--gs'); 
hold on;

%adding plots of Reference(ALL) (without CP)
plot(1:number_of_ranks,PSNR_vec_NoCp_ALL);
hold on;

legend('PSNR',['PSNR BEST= (',num2str(CP_Rank_Gaussian_1),' , ',num2str(maxPSNR_Gaussian_1_val_ALL),')'],...
   ['PSNR NoCp = ',num2str(PSNR_ALL_NoCP)]);

%printf(['the CP-rank is ' , num2str(CP_Rank) ])


%% Decomposition of the noised (gaussian 0.03) projection
% collect the data with many CP rank

TensorProj_Gaussian_3_CP_ALL = cell(1,number_of_ranks);
TensorProj_Clean_ALL = repmat({double(Tensor_Clean_Proj)},1,number_of_ranks);

% Starting the waitbar
counter = ConsoleProgressBar();
% Starting main iteration of MFISTA
counter.start(); counter.setText('Finding the best CP Decomposition Rank..');

for rank = 1:number_of_ranks
    TensorProj_Gaussian_3_CP_ALL{1,rank} = double(cp_als(Tensor_Gaussian_3_Proj,rank));
    % Collect data of the PSNR between the Tensor+Noise+CP and the clean Tensor
    PSNR_vec_Gaussian_3_ALL(rank) = psnr(TensorProj_Gaussian_3_CP_ALL{1,rank},TensorProj_Clean_ALL{1,rank});
    counter.setValue(rank/number_of_ranks);
end
counter.stop();

%find the max PSNR (value and index)
[maxPSNR_Gaussian_3_val_ALL , CP_Rank_Gaussian_3] = max(PSNR_vec_Gaussian_3_ALL(:));
%% Graph Section
%collect data for the comparing
%Calc the /PSNR/ of ALL proj
PSNR_ALL_NoCP = psnr(double(Tensor_Gaussian_3_Proj),double(Tensor_Clean_Proj));
% create /PSNR/ VECTOR (without CP) for the graph of ALL
PSNR_vec_NoCp_ALL = repmat(PSNR_ALL_NoCP,1,number_of_ranks);

% Create Graph of the Compare Between The Ground Truth and the Projections with several CP RANK.
figure(17);

%ALL projection
plot(1:number_of_ranks,PSNR_vec_Gaussian_3_ALL);
title('PSNR of the Proj after different decompositions');
xlim([0 number_of_ranks]);
xlabel('Rank');
ylabel('PSNR');
grid on;
hold on;

% mark the best RMSE/PSNR/SSIM  
plot(CP_Rank_Gaussian_3 , maxPSNR_Gaussian_3_val_ALL,'--gs'); 
hold on;

%adding plots of Reference(ALL) (without CP)
plot(1:number_of_ranks,PSNR_vec_NoCp_ALL);
hold on;

legend('PSNR',['PSNR BEST= (',num2str(CP_Rank_Gaussian_3),' , ',num2str(maxPSNR_Gaussian_3_val_ALL),')'],...
   ['PSNR NoCp = ',num2str(PSNR_ALL_NoCP)]);

%printf(['the CP-rank is ' , num2str(CP_Rank) ])
%% Decomposition of the noised (gaussian 0.05) projection
% collect the data with many CP rank

TensorProj_Gaussian_5_CP_ALL = cell(1,number_of_ranks);
TensorProj_Clean_ALL = repmat({double(Tensor_Clean_Proj)},1,number_of_ranks);

% Starting the waitbar
counter = ConsoleProgressBar();
% Starting main iteration of MFISTA
counter.start(); counter.setText('Finding the best CP Decomposition Rank..');

for rank = 1:number_of_ranks
    TensorProj_Gaussian_5_CP_ALL{1,rank} = double(cp_als(Tensor_Gaussian_5_Proj,rank));
    % Collect data of the PSNR between the Tensor+Noise+CP and the clean Tensor
    PSNR_vec_Gaussian_5_ALL(rank) = psnr(TensorProj_Gaussian_5_CP_ALL{1,rank},TensorProj_Clean_ALL{1,rank});
    counter.setValue(rank/number_of_ranks);
end
counter.stop();
%find the max PSNR (value and index)
[maxPSNR_Gaussian_5_val_ALL , CP_Rank_Gaussian_5] = max(PSNR_vec_Gaussian_5_ALL(:));

%% Graph Section
%collect data for the comparing
%Calc the /PSNR/ of ALL proj
PSNR_ALL_NoCP = psnr(double(Tensor_Gaussian_5_Proj),double(Tensor_Clean_Proj));
% create /PSNR/ VECTOR (without CP) for the graph of ALL
PSNR_vec_NoCp_ALL = repmat(PSNR_ALL_NoCP,1,number_of_ranks);

% Create Graph of the Compare Between The Ground Truth and the Projections with several CP RANK.
figure(17);

%ALL projection
plot(1:number_of_ranks,PSNR_vec_Gaussian_5_ALL);
title('PSNR of the Proj after different decompositions');
xlim([0 number_of_ranks]);
xlabel('Rank');
ylabel('PSNR');
grid on;
hold on;

% mark the best RMSE/PSNR/SSIM  
plot(CP_Rank_Gaussian_5 , maxPSNR_Gaussian_5_val_ALL,'--gs'); 
hold on;

%adding plots of Reference(ALL) (without CP)
plot(1:number_of_ranks,PSNR_vec_NoCp_ALL);
hold on;

legend('PSNR',['PSNR BEST= (',num2str(CP_Rank_Gaussian_5),' , ',num2str(maxPSNR_Gaussian_5_val_ALL),')'],...
   ['PSNR NoCp = ',num2str(PSNR_ALL_NoCP)]);

%printf(['the CP-rank is ' , num2str(CP_Rank) ])