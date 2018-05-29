%% Spectral CT
% Project by Gil Ben Ari and Rotem Zamir
 clear
 close all;
 clc;
% to run script you need to add all the CT folder to the path.


%% Spectra Initialization
%Proj Tensors Load
%load('P:\Files\TensorProj_MatArray.mat');
%load('P:\Files\TensorProjNoiseCP_MatArray');

% Number of spectra involved in each scan
SpectralParams.K = 64;

% Amount of photons emitted, proportional to dosage
SpectralParams.N0 = 10^20;

% The tube voltage, if we are on "switching" it should have length(K) elements
% else (thresholding) just the center frequency of the single emitted spectrum
% The Vtube range is [30 140] 
SpectralParams.Vtube = 130;
%SpectralParams.Vtube = 100;

% The bin right limits markers, the first bin starts at 0   
SpectralParams.Threshold = 67:130; %for 64 enegies 
%SpectralParams.Threshold = [50 60 70 80 90 100];

% The analog energy interpolation factor for the integral approximation
SpectralParams.InterpolationRes = 1000;

% The ground truth phantom side resolution (NxN phantom, N=PhantomRes)
PhantomRes = 64;

%sigma of gaussian noise
sigma = 0.03;

%number of rank for the CP simulation
number_of_ranks = 150;

SpectralParams.ToOptimizeThresholds = 1; 

%% choosing materials
% materials for MusFromMaterials func (will be inside BodyMaterials):
% adipose, blood, bone_compact, bone_cortical, brain, lung, muscle_skeletal,
% muscle_striated, skin, soft_tissue, water

% materials for XrayMu func will be string with Chemical formula or simple
% atomic notation

I = 3; %Number of materials
ChosenMaterials = { 'I','Ca', 'soft_tissue'}; %choose the materials you wish for the basis
% FIXME: A normalizing vector, trying to equalize energies between materials
base = [0.125 1 20]./255;
%base = [1 1 1];

%% Loading phantom
% chosen phantom can be either: 'circles', 'more_circles', 'paint',
% 'spectral_shepp', 'spectral_zubal'.
chosen_phantom = 'circles';
RGBphantom = LoadPhantom(PhantomRes,chosen_phantom);
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
%Spectra = GetSpectra('Thresholding', SpectralParams);
Spectra = GetSpectra('Thresholding', SpectralParams, ChosenMaterials, base); 

%% Creating the scan of the object (simulating a spectral scan)
[ Spectra, projections, RadonMaterials, ChosenMaterials ] = scan( Spectra, phantom, ChosenMaterials );

%% Defining the Tensors from the scan cell data
ProjectionTensor = zeros([size(projections{1}),length(projections)]);
for i=1:length(projections)
    ProjectionTensor(:,:,i) = projections{i};
end

%% Defining the Projection Tensor With the Tensor-ToolBox and display with differents ranks.

TensorProj=tensor(ProjectionTensor); %% define the tensor with the toolbox

%% Defining the Projection Tensor With Gussian Noise
ProjectionTensorWithNoise = zeros([size(projections{1}),length(projections)]);
for i=1:SpectralParams.K
    % GaussianNoise
    [ProjectionTensorWithNoise(:,:,i),~] = GaussianNoise( ProjectionTensor(:,:,i),sigma );
    % simple noise
    %ProjectionTensorWithNoise(:,:,i) = ProjectionTensor(:,:,i) + sigma.*range(ProjectionTensor(:,:,i)).*randn(size(ProjectionTensor(:,:,i)));
end

%% Check and Display the RMSE of each Projection with noise and the RMSE of each Projection with noise and CP.
TensorProjNoise=tensor(ProjectionTensorWithNoise);

%% Collect Data for Compare Between The Ground Truth and the Projections with several CP RANK.
% Extremly slow script (2.5 hours of running). I save the data of the script after I run this
% script once at P:\Files\.

% collect the data of the full image (all the projections) with many CP rank
TensorProjNoiseCP_ALL = cell(1,number_of_ranks);
TensorProj_ALL = repmat({double(TensorProj)},1,number_of_ranks);

for rank = 1:number_of_ranks
    TensorProjNoiseCP_ALL{1,rank} = double(cp_als(TensorProjNoise,rank));
end


%% Compare Between The Ground Truth and the Projections with several CP RANK.

% Compare full image - Ground Truth and all the Projections with several CP RANK.
for rank = 1:number_of_ranks
    RMSE_vec_ALL(rank) = rmse(TensorProjNoiseCP_ALL{1,rank},TensorProj_ALL{1,rank});
    PSNR_vec_ALL(rank) = psnr(TensorProjNoiseCP_ALL{1,rank},TensorProj_ALL{1,rank});
    [SSIM_vec_ALL(rank), ~] = ssim(TensorProjNoiseCP_ALL{1,rank},TensorProj_ALL{1,rank});
end

%Calc the RMSE/PSNR/SSIM of ALL proj
RMSE_ALL_NoCP = rmse(double(TensorProjNoise),double(TensorProj));
PSNR_ALL_NoCP = psnr(double(TensorProjNoise),double(TensorProj));
[SSIM_ALL_NoCP, ~] = ssim(double(TensorProjNoise),double(TensorProj));

% create RMSE/PSNR/SSIM VECTOR (without CP) for the graph of ALL
RMSE_vec_NoCp_ALL = repmat(RMSE_ALL_NoCP,1,number_of_ranks);
PSNR_vec_NoCp_ALL = repmat(PSNR_ALL_NoCP,1,number_of_ranks);
SSIM_vec_NoCp_ALL = repmat(SSIM_ALL_NoCP,1,number_of_ranks);

%% finding the max of the PSNR and the min of the RMSE for ALL projection
[minRMSEval_ALL,minRMSEindex_ALL] = min(RMSE_vec_ALL(:));
[maxPSNRval_ALL,maxPSNRindex_ALL] = max(PSNR_vec_ALL(:));
[maxSSIMval_ALL,maxSSIMindex_ALL] = max(SSIM_vec_ALL(:));

%% Create Graph of the Compare Between The Ground Truth and the each Projection with several CP RANK.
RANK = 1:number_of_ranks;
figure(9);

%ALL projection
%subplot(4,1,4);
plot(1:number_of_ranks,RMSE_vec_ALL,'o');
title(['Compare Between The Ground Truth and ', '\bf ALL ', '\rm of the Projection with several CP RANK.']);
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
legend('RMSE','PSNR','SSIM',['RMSE BEST= (',num2str(minRMSEindex_ALL),' , ',num2str(minRMSEval_ALL),')'],...
    ['PSNR BEST= (',num2str(maxPSNRindex_ALL),' , ',num2str(maxPSNRval_ALL),')'],...
    ['SSIM BEST= (',num2str(maxSSIMindex_ALL),' , ',num2str(maxSSIMval_ALL),')'],...
    ['RMSE NoCp = ',num2str(RMSE_ALL_NoCP)],['PSNR NoCp = ',num2str(PSNR_ALL_NoCP)],['SSIM NoCp = ',num2str(SSIM_ALL_NoCP)]);

% Label: ALL 16_8.mat


%% Reconstrucion Analyze 

% Convert the Tensors to cell of matrixes
% Defining the projections cell from the projections Tensor (before CP).
Projections_Noise =  TensorToMatrixInCell(TensorProjNoise,6);  
% Defining the projections cell from the projections Tensor (after CP) for all the ranks.
Projections_Noise_CP_RankCell = cell(1,number_of_ranks);
for i = 1:number_of_ranks
    % TensorProjNoiseCP_ALL is the Tensor cell (define in previouses sections). 
    % in each{1,rank} there is the CP with specific rank.
    Projections_Noise_CP_RankCell{i} =  TensorToMatrixInCell(TensorProjNoiseCP_ALL{1,i},6);  
end


%% Collect data of Reconstruction of the noised projection with CP

estimate_phantom = cell(SpectralParams.K , number_of_ranks);

for i = 1:number_of_ranks
    temp_Projections_Noise_Rank = Projections_Noise_CP_RankCell{i}; %temp_Projections_Noise_Rank is the (Kx1)Cell of specific rank (k projections)
    for j = 1:SpectralParams.K 
        % estimate_phantom{j,i} will hold the j Projection with noise, after
        % CP with rank i
        % conjgrad is function from the SpectralSimulator/Spectral_FISTA_recounstruction (line 110)  
        % of Rotem%Gill (function of Shahar). its get matrix of Projections (in our case temp_Projections_Noise_Rank{j})
        % and return estimate phantom, the 'reconstructed image'. I think...
        estimate_phantom{j,i} = conjgrad(@(x) Rpp_T(M(Rpp(x))) , Rpp_T(M(temp_Projections_Noise_Rank{j})),zeros(128),50,1e-60);
    end
end

%% Collect data of Reconstruction of the noised projection without CP

estimate_phantom_withoutCP = cell(SpectralParams.K , 1);

for j = 1:SpectralParams.K 
    % estimate_phantom_withoutCP{j,1} will hold the j Projection with noise
    estimate_phantom_withoutCP{j,1} = conjgrad(@(x) Rpp_T(M(Rpp(x))) , Rpp_T(M(Projections_Noise{j})),zeros(128),50,1e-60);
end

%% Collect data of Reconstruction of the clean projection (without CP)

estimate_phantom_clean = cell(SpectralParams.K , 1);

for j = 1:SpectralParams.K 
    % estimate_phantom_withoutCP{j,1} will hold the j Projection with noise
    estimate_phantom_clean{j,1} = conjgrad(@(x) Rpp_T(M(Rpp(x))) , Rpp_T(M(projections{j})),zeros(128),50,1e-60);
end


%% Collect Data of the RMSE/PSNR/SSIM of the estimate phantom

% Compare Between The Ground Truth and the estimate phantom with several CP RANK for each projection.

proj_num=1;
RMSE_Reconstructed_mat=0;
PSNR_Reconstructed_mat=0;
SSIM_Reconstructed_mat=0; 
RMSE_Reconstructed_NoCp_vec = cell(SpectralParams.K , 1);
PSNR_Reconstructed_NoCp_vec = cell(SpectralParams.K , 1);
SSIM_Reconstructed_NoCp_vec = cell(SpectralParams.K , 1);


for proj_num = 1:SpectralParams.K
    % Compare between each rank of each reconstructed projection with noise
    % with the clean reconstructed projection
    for rank = 1:number_of_ranks
        y = estimate_phantom{proj_num , rank};
        y_ref = estimate_phantom_clean{proj_num,1};
        RMSE_Reconstructed_mat(rank,proj_num) = rmse(y,y_ref);
        PSNR_Reconstructed_mat(rank,proj_num) = psnr(y,y_ref);
        [SSIM_Reconstructed_mat(rank,proj_num), ~] = ssim(y,y_ref);
    end

    %Compare Between The Ground Truth and the Reconstructed Projections (with the noise) without CP.
    y = estimate_phantom_withoutCP{proj_num , 1};
    y_ref = estimate_phantom_clean{proj_num,1};
    RMSE_Reconstructed_NoCp(proj_num) = rmse(y,y_ref);
    PSNR_Reconstructed_NoCp(proj_num) = psnr(y,y_ref);
    [SSIM_Reconstructed_NoCp(proj_num), ~] = ssim(y,y_ref);
    
    % create RMSE/PSNR/SSIM VECTOR (without CP) for the graph of each proj
    RMSE_Reconstructed_NoCp_vec{proj_num,1} = repmat(RMSE_Reconstructed_NoCp(proj_num),1,number_of_ranks);
    PSNR_Reconstructed_NoCp_vec{proj_num,1} = repmat(PSNR_Reconstructed_NoCp(proj_num),1,number_of_ranks);
    SSIM_Reconstructed_NoCp_vec{proj_num,1} = repmat(SSIM_Reconstructed_NoCp(proj_num),1,number_of_ranks);
end


%% finding the max of the PSNR and the min of the RMSE for each Reconstructed projection
for proj_num = 1:SpectralParams.K
    [minRMSEval(proj_num),minRMSEindex(proj_num)] = min(RMSE_Reconstructed_mat(:,proj_num));
    [maxPSNRval(proj_num),maxPSNRindex(proj_num)] = max(PSNR_Reconstructed_mat(:,proj_num));
    [maxSSIMval(proj_num),maxSSIMindex(proj_num)] = max(SSIM_Reconstructed_mat(:,proj_num));
end

% finding the max of the PSNR and the min of the RMSE for ALL projection
[minRMSEval_ALL,minRMSEindex_ALL] = min(RMSE_vec_ALL(:));
[maxPSNRval_ALL,maxPSNRindex_ALL] = max(PSNR_vec_ALL(:));
[maxSSIMval_ALL,maxSSIMindex_ALL] = max(SSIM_vec_ALL(:));




%% Create Graph of the Compare Between The Ground Truth and the each Reconstructed Projection with several CP RANK.
RANK = 1:number_of_ranks;
figure(14);
for i = 1:SpectralParams.K
    subplot(3,2,i);
    %i=1;
    plot(RANK,RMSE_Reconstructed_mat(:,i),'o');
    title(['Compare Between The Ground Truth and the ', num2str(i) , ' Reconstructed Projection with several CP RANK.']);
    xlim([0 100]);
    grid on;
    hold on;
    plot(RANK,PSNR_Reconstructed_mat(:,i),'o');
    hold on;
    plot(RANK,SSIM_Reconstructed_mat(:,i),'o');
    % mark the best RMSE/PSNR/SSIM
    plot(minRMSEindex(i) , minRMSEval(i),'--gs');
    hold on;
    plot(maxPSNRindex(i) , maxPSNRval(i),'--gs'); 
    hold on;
    plot(maxSSIMindex(i) , maxSSIMval(i),'--gs'); 
    hold on;
    % adding plots of Reference (without CP)
    plot(1:number_of_ranks,RMSE_Reconstructed_NoCp_vec{i,1});
    hold on;
    plot(1:number_of_ranks,PSNR_Reconstructed_NoCp_vec{i,1});
    hold on;
    plot(1:number_of_ranks,SSIM_Reconstructed_NoCp_vec{i,1});
    hold on;
    legend('RMSE','PSNR','SSIM',['RMSE_BEST= (',num2str(minRMSEindex(i)),' , ',num2str(minRMSEval(i)),')'],['PSNR_BEST= (',num2str(maxPSNRindex(i)),' , ',num2str(maxPSNRval(i)),')'],['SSIM_BEST= (',num2str(maxSSIMindex(i)),' , ',num2str(maxSSIMval(i)),')'],...
        ['RMSE_NoCp = ',num2str(RMSE_Reconstructed_NoCp_vec{i,1}(1))],['PSNR_NoCp = ',num2str(PSNR_Reconstructed_NoCp_vec{i,1}(1))],['SSIM_NoCp = ',num2str(SSIM_Reconstructed_NoCp_vec{i,1}(1))]);
end

