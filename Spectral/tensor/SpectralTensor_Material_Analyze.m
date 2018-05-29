%% Spectral CT
% Project by Gil Ben Ari and Rotem Zamir
 clear
 close all;
 clc;
% to run script you need to add all the CT folder to the path.


%% Spectra Initialization

% Number of spectra involved in each scan
SpectralParams.K = 6;

% Amount of photons emitted, proportional to dosage
SpectralParams.N0 = 10^20;

% The tube voltage, if we are on "switching" it should have length(K) elements
% else (thresholding) just the center frequency of the single emitted spectrum
% The Vtube range is [30 140]
SpectralParams.Vtube = 100;

% The bin right limits markers, the first bin starts at 0
SpectralParams.Threshold = [ 50,60,70,80,90,100];

% The analog energy interpolation factor for the integral approximation
SpectralParams.InterpolationRes = 1000;

% The ground truth phantom side resolution (NxN phantom, N=PhantomRes)
PhantomRes = 64;

%sigma of gaussian noise
sigma = 0.03;

%number of rank for the CP simulation
number_of_ranks = 100;

%% choosing materials
% materials for MusFromMaterials func (will be inside BodyMaterials):
% adipose, blood, bone_compact, bone_cortical, brain, lung, muscle_skeletal,
% muscle_striated, skin, soft_tissue, water

% materials for XrayMu func will be string with Chemical formula or simple
% atomic notation

I = 6; %Number of materials
ChosenMaterials = { 'I','Ca', 'soft_tissue','blood','bone_compact','muscle_skeletal'}; %choose the materials you wish for the basis
% FIXME: A normalizing vector, trying to equalize energies between materials
% base = [0.125 1 20]./255;
base = [0.125 1 20 10 5 10 ]./255;

%% Loading phantom
% chosen phantom can be either: 'circles', 'more_circles', 'paint',
% 'spectral_shepp', 'spectral_zubal'.
chosen_phantom = 'circles';
RGBphantom1 = LoadPhantom(PhantomRes,chosen_phantom);
RGBphantom2 = LoadPhantom(PhantomRes,chosen_phantom);
RGBphantom = cat(3,RGBphantom1,RGBphantom2);
phantom = cell(I,1);

for ii=1:I
    phantom{ii} = base(ii)*RGBphantom(:,:,ii);
end


%% Creating spectra
%Spectra:
%scantype defines the type of the scan for the simulation. It can be:
%'switching' - used mainly for DECT, where the tube's voltage changes (usually 2 different
% voltages). This option simulates the needed spectra.
%'Thresholding' - used for spectral scan, where the detectors are energy
%sensitive. This option simulates the thresholded bins' spectra.
%'ExtSpectra' - used when there is an external measured spectra, and want to
% simulate the corresponding attenuation

% calling GetSpectra or loading external spectra data
Spectra = GetSpectra('Thresholding', SpectralParams);


%% Creating the scan of the object (simulating a spectral scan)
[ Spectra, projections, RadonMaterials, ChosenMaterials ] = scan( Spectra, phantom, ChosenMaterials );


%% Defining the Tensors from the scan cell data
MaterialTensor   = zeros([size(RadonMaterials{1}),length(RadonMaterials)]);
for i=1:length(RadonMaterials)
    MaterialTensor(:,:,i) = RadonMaterials{i};
end


%% Defining the Material Tensor With the Tensor-ToolBox and display with differents ranks.
TensorMaterial=tensor(MaterialTensor); %% define the tensor with the toolbox

%% Defining the Material Tensor With Gussian Noise

MaterialTensorWithNoise = zeros([size(RadonMaterials{1}),length(RadonMaterials)]);
for i=1:I
    [MaterialTensorWithNoise(:,:,i),~] = GaussianNoise( MaterialTensor(:,:,i),sigma );
end

%% Check and Display the RMSE of each Material with noise and the RMSE of each Material with noise and CP.
TensorMaterialNoise=tensor(MaterialTensorWithNoise);

%% Collect Data for Compare Between The Ground Truth and the Projections with several CP RANK.

% collect the data of the full image (all the projections) with many CP rank
TensorMaterialNoiseCP_ALL = cell(1,number_of_ranks);
TensorMaterial_ALL = repmat({double(TensorMaterial)},1,number_of_ranks);

for rank = 1:number_of_ranks
    TensorMaterialNoiseCP_ALL{1,rank} = double(cp_als(TensorMaterialNoise,rank));
end


%% Compare Between The Ground Truth and the Projections with several CP RANK.

% Compare full image - Ground Truth and all the Projections with several CP RANK.
for rank = 1:number_of_ranks
    RMSE_vec_ALL(rank) = rmse(TensorMaterialNoiseCP_ALL{1,rank},TensorMaterial_ALL{1,rank});
    PSNR_vec_ALL(rank) = psnr(TensorMaterialNoiseCP_ALL{1,rank},TensorMaterial_ALL{1,rank});
    [SSIM_vec_ALL(rank), ~] = ssim(TensorMaterialNoiseCP_ALL{1,rank},TensorMaterial_ALL{1,rank});
end

%Calc the RMSE/PSNR/SSIM of ALL proj
RMSE_ALL_NoCP = rmse(double(TensorMaterialNoise),double(TensorMaterial));
PSNR_ALL_NoCP = psnr(double(TensorMaterialNoise),double(TensorMaterial));
[SSIM_ALL_NoCP, ~] = ssim(double(TensorMaterialNoise),double(TensorMaterial));

% create RMSE/PSNR/SSIM VECTOR (without CP) for the graph of ALL
RMSE_vec_NoCp_ALL = repmat(RMSE_ALL_NoCP,1,number_of_ranks);
PSNR_vec_NoCp_ALL = repmat(PSNR_ALL_NoCP,1,number_of_ranks);
SSIM_vec_NoCp_ALL = repmat(SSIM_ALL_NoCP,1,number_of_ranks);

%% finding the max of the PSNR and the min of the RMSE for ALL projection
[minRMSEval_ALL,minRMSEindex_ALL] = min(RMSE_vec_ALL(:));
[maxPSNRval_ALL,maxPSNRindex_ALL] = max(PSNR_vec_ALL(:));
[maxSSIMval_ALL,maxSSIMindex_ALL] = max(SSIM_vec_ALL(:));

%% Create Graph of the Compare Between The Ground Truth and the each Material with several CP RANK.
RANK = 1:number_of_ranks;
figure(9);

%ALL Material
%subplot(4,1,4);
plot(1:number_of_ranks,RMSE_vec_ALL,'o');
title(['Compare Between The Ground Truth and ', '\bf ALL ', '\rm of the Material with several CP RANK.']);
xlim([0 number_of_ranks]);
grid on;
hold on;
plot(1:number_of_ranks,PSNR_vec_ALL,'o');
hold on;
plot(1:number_of_ranks,SSIM_vec_ALL,'o');
% mark the best RMSE/PSNR/SSIM  
plot(minRMSEindex_ALL , minRMSEval_ALL,'--gs');
hold on;
plot(maxPSNRindex_ALL , maxPSNRval_ALL,'--gs'); 
hold on;
plot(maxSSIMindex_ALL , maxSSIMval_ALL,'--gs'); 
hold on;
%adding plots of Reference(ALL) (without CP)
plot(1:number_of_ranks,RMSE_vec_NoCp_ALL);
hold on;
plot(1:number_of_ranks,PSNR_vec_NoCp_ALL);
hold on;
plot(1:number_of_ranks,SSIM_vec_NoCp_ALL);
hold on;
legend('RMSE','PSNR','SSIM',['RMSE_BEST= (',num2str(minRMSEindex_ALL),' , ',num2str(minRMSEval_ALL),')'],...
    ['PSNR_BEST= (',num2str(maxPSNRindex_ALL),' , ',num2str(maxPSNRval_ALL),')'],...
    ['SSIM_BEST= (',num2str(maxSSIMindex_ALL),' , ',num2str(maxSSIMval_ALL),')'],...
    ['RMSE_NoCp = ',num2str(RMSE_ALL_NoCP)],['PSNR_NoCp = ',num2str(PSNR_ALL_NoCP)],['SSIM_NoCp = ',num2str(SSIM_ALL_NoCP)]);