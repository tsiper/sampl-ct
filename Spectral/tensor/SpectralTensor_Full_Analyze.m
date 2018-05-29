%% Spectral CT
% Project by Gil Ben Ari and Rotem Zamir
 clear
 close all;
 clc;
% to run script you need to add all the CT folder to the path.


%% Spectra Initialization
%Proj Tensors Load
load('P:\Files\TensorProj_MatArray.mat');
load('P:\Files\TensorProjNoiseCP_MatArray');

% Number of spectra involved in each scan
SpectralParams.K = 3;

% Amount of photons emitted, proportional to dosage
SpectralParams.N0 = 10^20;

% The tube voltage, if we are on "switching" it should have length(K) elements
% else (thresholding) just the center frequency of the single emitted spectrum
% The Vtube range is [30 140]
SpectralParams.Vtube = [100];

% The bin right limits markers, the first bin starts at 0
SpectralParams.Threshold = [ 50, 75, 100];

% The analog energy interpolation factor for the integral approximation
SpectralParams.InterpolationRes = 1000;

% The ground truth phantom side resolution (NxN phantom, N=PhantomRes)
PhantomRes = 64;

%% choosing materials
% materials for MusFromMaterials func (will be inside BodyMaterials):
% adipose, blood, bone_compact, bone_cortical, brain, lung, muscle_skeletal,
% muscle_striated, skin, soft_tissue, water

% materials for XrayMu func will be string with Chemical formula or simple
% atomic notation

I = 3; %Number of materials
ChosenMaterials = { 'I','Ca', 'soft_tissue'}; %choose the materials you wish for the basis
% FIXME: A normalizing vector, trying to equalize energies between materials
% base = [0.125 1 20]./255;
base = [1 1 1];

%% Loading phantom
% chosen phantom can be either: 'circles', 'more_circles', 'paint',
% 'spectral_shepp', 'spectral_zubal'.
chosen_phantom = 'circles';
RGBphantom = LoadPhantom(PhantomRes,chosen_phantom);
phantom = cell(I,1);

for ii=1:I
    phantom{ii} = base(ii)*RGBphantom(:,:,ii);
end

%% Thresholds optimization by Alvarez


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

%%

%NEED TO DO: Normalize phantom levels to actual range of thicknesses (graylevel = 255);

%% Creating the scan of the object (simulating a spectral scan)
[ Spectra, projections, RadonMaterials, ChosenMaterials ] = scan( Spectra, phantom, ChosenMaterials );


%% Defining the Tensors from the scan cell data
MaterialTensor   = zeros([size(RadonMaterials{1}),length(RadonMaterials)]);
ProjectionTensor = zeros([size(projections{1}),length(projections)]);
for i=1:length(RadonMaterials)
    MaterialTensor(:,:,i) = RadonMaterials{i};
end
for i=1:length(projections)
    ProjectionTensor(:,:,i) = projections{i};
end
figure;
subplot(121);imagesc(ProjectionTensor(:,:,1:3));
title('Projection Tensor');
subplot(122);imagesc(MaterialTensor);
title('Material Tensor');

%% Defining the Projection Tensor With the Tensor-ToolBox and display with differents ranks.

TensorProj=tensor(ProjectionTensor); %% define the tensor with the toolbox
figure(1); 
j=2;
subplot(5,5,1); imagesc(ProjectionTensor(:,:,1:3));title('the Projection Tensor'); %subplot(5,5,2); imagesc(double(TensorProj));title('the Tensor projection');

for i=1:10
    subplot(5,5,j);imagesc(double(cp_als(TensorProj,i))); title(['CP with Rank = ', num2str(i) ,' of Tensor Projection']);
    j=j+1; 
end

for i=15:5:60
    subplot(5,5,j);imagesc(double(cp_als(TensorProj,i))); title(['CP with Rank = ', num2str(i) ,' of Tensor Projection']);
    j=j+1; 
end

for i = 70:100:300
    subplot(5,5,j);imagesc(double(cp_als(TensorProj,i))); title(['CP with Rank = ', num2str(i) ,' of Tensor Projection']);
    j=j+1; 
end
i=1000;
subplot(5,5,j);imagesc(double(cp_als(TensorProj,i))); title(['CP with Rank = ', num2str(i) ,' of Tensor Projection']);


%% Display each projection in the Tensor Projection
figure(2); 
rank = 80;
j=1;
for i=1:3
    subplot(3,2,j);imagesc(double(cp_als(TensorProj(:,:,i),rank))); title(['CP with Rank = ', num2str(rank) ,' of #', num2str(i),' Projection']);
    j=j+1;
    subplot(3,2,j); imagesc(ProjectionTensor(:,:,i)); title(['The of #', num2str(i),' Projection']);
    j=j+1;
end

%% Defining the Projection Tensor With Gussian Noise
ProjectionTensorWithNoise = zeros([size(projections{1}),length(projections)]);
sigma = 0.03;
for i=1:3
    [ProjectionTensorWithNoise(:,:,i),~] = GaussianNoise( ProjectionTensor(:,:,i),sigma );
end

%% Display the Projection Tensor With Gussian Noise for each projection
figure(3);
subplot(221);imagesc(ProjectionTensorWithNoise(:,:,1:3)); title('Projection Tensor With Noise'); 
for i=1:3
    subplot(2,2,(i+1));imagesc(ProjectionTensorWithNoise(:,:,i)); title(['The #',num2str(i),' Projection With Noise']); 
end

%% Check and Display the RMSE of each Projection with noise and the RMSE of each Projection with noise and CP.
TensorProjNoise=tensor(ProjectionTensorWithNoise);
figure(4);
rank = 23;
subplot(4,4,1); imagesc(double(TensorProjNoise)); title('Tensor Projection with noise');
RMSE = rmse(double(TensorProjNoise),double(TensorProj));
subplot(4,4,2); imagesc(double(TensorProjNoise)-double(TensorProj)); title(['The MSE of Projection Tensor with noise = ', num2str(RMSE) ]);
subplot(4,4,3); imagesc(double(cp_als(TensorProjNoise,rank))); title(['CP with Rank = ', num2str(rank) ,' of Tensor Projection with noise']);
RMSE = rmse(double(cp_als(TensorProjNoise,rank)),double(TensorProj));
subplot(4,4,4); imagesc(double(cp_als(TensorProjNoise,rank))- double(TensorProj)); title(['MSE of CP with Rank ', num2str(rank) ,' of Projection with noise = ',num2str(RMSE)]);

j=5;
for i = 1:3
    subplot(4,4,j); imagesc(double(TensorProjNoise(:,:,i))); title(['The #',num2str(i),' Projection with noise']);
    j=j+1;
    RMSE = rmse(double(TensorProjNoise(:,:,i)),double(TensorProj(:,:,i)));
    subplot(4,4,j); imagesc(double(TensorProjNoise(:,:,i)) - double(TensorProj(:,:,i))); title(['The MSE of #',num2str(i),' proj with noise = ', num2str(RMSE) ]);
    j=j+1;
    
    subplot(4,4,j); imagesc(double(cp_als(TensorProjNoise(:,:,i),rank))); title(['CP with Rank = ', num2str(rank) ,' of #', num2str(i),' Projection with noise']);
    j=j+1;
    RMSE = rmse(double(cp_als(TensorProjNoise(:,:,i),rank)),double(TensorProj(:,:,i)));
    subplot(4,4,j); imagesc(double(cp_als(TensorProjNoise(:,:,i),rank)) - double(TensorProj(:,:,i))); title(['MSE of CP Rank = ', num2str(rank) ' of #',num2str(i),' proj with noise = ', num2str(RMSE) ]);
    j=j+1;
end

%% Defining the Material Tensor With the Tensor-ToolBox and display with differents ranks.

TensorMaterial=tensor(MaterialTensor); %% define the tensor with the toolbox
figure(5); 
j=2;
subplot(5,5,1); imagesc(double(TensorMaterial(:,:,1:3)));title('the Material Tensor'); %subplot(5,5,2); imagesc(double(TensorMaterial));title('the Tensor Material');

for i=1:10
    subplot(5,5,j);imagesc(double(cp_als(TensorMaterial,i))); title(['CP with Rank = ', num2str(i) ,' of Material Tensor']);
    j=j+1; 
end

for i=15:5:60
    subplot(5,5,j);imagesc(double(cp_als(TensorMaterial,i))); title(['CP with Rank = ', num2str(i) ,' of Material Tensor']);
    j=j+1; 
end

for i = 70:100:300
    subplot(5,5,j);imagesc(double(cp_als(TensorMaterial,i))); title(['CP with Rank = ', num2str(i) ,' of Material Tensor']);
    j=j+1; 
end
i=1000;
subplot(5,5,j);imagesc(double(cp_als(TensorMaterial,i))); title(['CP with Rank = ', num2str(i) ,' of Material Tensor']);

%% Display each Material in the Material Tensor 
figure(6); 
rank = 80;
j=1;
for i=1:3
    subplot(3,2,j);imagesc(double(cp_als(TensorMaterial(:,:,i),rank))); title(['CP with Rank = ', num2str(rank) ,' of #', num2str(i),' Material']);
    j=j+1;
    subplot(3,2,j); imagesc(MaterialTensor(:,:,i)); title(['The of #', num2str(i),' Material']);
    j=j+1;
end

%% Defining the Material Tensor With Gussian Noise

MaterialTensorWithNoise = zeros([size(RadonMaterials{1}),length(RadonMaterials)]);
sigma = 0.03;
for i=1:3
    [MaterialTensorWithNoise(:,:,i),~] = GaussianNoise( MaterialTensor(:,:,i),sigma );
end
 
%% Display the Material Tensor With Gussian Noise for each Material
figure(7);
subplot(221);imagesc(MaterialTensorWithNoise(:,:,1:3)); title('Material Tensor With Noise'); 
for i=1:3
    subplot(2,2,(i+1));imagesc(MaterialTensorWithNoise(:,:,i)); title(['The #',num2str(i),' Material With Noise']); 
end

%% Check and Display the RMSE of each Material with noise and the RMSE of each Material with noise and CP.
TensorMaterialNoise=tensor(MaterialTensorWithNoise);
figure(8);
rank = 80;
subplot(4,4,1); imagesc(double(TensorMaterialNoise)); title('Tensor Material with noise');
RMSE = rmse(double(TensorMaterialNoise),double(TensorMaterial));
subplot(4,4,2); imagesc(double(TensorMaterialNoise)-double(TensorMaterial)); title(['The MSE of Material Tensor with noise = ', num2str(RMSE) ]);
subplot(4,4,3); imagesc(double(cp_als(TensorMaterialNoise,rank))); title(['CP with Rank = ', num2str(rank) ,' of Tensor Material with noise']);
RMSE = rmse(double(cp_als(TensorMaterialNoise,rank)),double(TensorMaterial));
subplot(4,4,4); imagesc(double(cp_als(TensorMaterialNoise,rank))- double(TensorMaterial)); title(['MSE of CP with Rank ', num2str(rank) ,' of Material with noise = ',num2str(RMSE)]);

j=5;
for i = 1:3
    subplot(4,4,j); imagesc(double(TensorMaterialNoise(:,:,i))); title(['The #',num2str(i),' Material with noise']);
    j=j+1;
    RMSE = rmse(double(TensorMaterialNoise(:,:,i)),double(TensorMaterial(:,:,i)));
    subplot(4,4,j); imagesc(double(TensorMaterialNoise(:,:,i)) - double(TensorMaterial(:,:,i))); title(['The MSE of #',num2str(i),' Material with noise = ', num2str(RMSE) ]);
    j=j+1;
    
    subplot(4,4,j); imagesc(double(cp_als(TensorMaterialNoise(:,:,i),rank))); title(['CP with Rank = ', num2str(rank) ,' of #', num2str(i),' Material with noise']);
    j=j+1;
    RMSE = rmse(double(cp_als(TensorMaterialNoise(:,:,i),rank)),double(TensorMaterial(:,:,i)));
    subplot(4,4,j); imagesc(double(cp_als(TensorMaterialNoise(:,:,i),rank)) - double(TensorMaterial(:,:,i))); title(['MSE of CP Rank = ', num2str(rank) ' of #',num2str(i),' Material with noise = ', num2str(RMSE) ]);
    j=j+1;
end


%% Collect Data for Compare Between The Ground Truth and the Projections with several CP RANK.
% Extremly slow script (2.5 hours of running). I save the data of the script after I run this
% script once at P:\Files\.

number_of_projections = 3;
number_of_ranks = 350;
TensorProjNoiseCP_MatArray = cell(number_of_ranks , number_of_projections);
TensorProj_MatArray = cell(number_of_ranks , number_of_projections);
for proj_num = 1:number_of_projections
i=1;
    for rank = 1:number_of_ranks
        TensorProjNoiseCP_MatArray{i,proj_num} = double(cp_als(TensorProjNoise(:,:,proj_num),rank));
        %there is option to improve the function. to get outside from the loop the next line. 
        TensorProj_MatArray{i,proj_num} = double(TensorProj(:,:,proj_num));
        i=i+1;
    end
end

% collect the data of the full image (all the projections) with many CP
% rank
TensorProjNoiseCP_ALL = cell(1,50);
TensorProj_ALL = repmat({double(TensorProj)},1,50);

for rank = 1:50
    TensorProjNoiseCP_ALL{1,rank} = double(cp_als(TensorProjNoise,rank));
end



% load('P:\Files\TensorProj_MatArray.mat');
% load('P:\Files\TensorProjNoiseCP_MatArray');

%% Compare Between The Ground Truth and the Projections with several CP RANK.
proj_num=1;
RMSE_mat=0;
PSNR_mat=0;
SSIM_mat=0; 
RMSE_NoCp_vec = cell(number_of_projections , 1);
PSNR_NoCp_vec = cell(number_of_projections , 1);
SSIM_NoCp_vec = cell(number_of_projections , 1);


for proj_num = 1:number_of_projections
    for rank = 1:number_of_ranks
        y = TensorProjNoiseCP_MatArray{rank,proj_num};
        y_ref = TensorProj_MatArray{rank,proj_num};
        RMSE_mat(rank,proj_num) = rmse(y,y_ref);
        PSNR_mat(rank,proj_num) = psnr(y,y_ref);
        [SSIM_mat(rank,proj_num), ~] = ssim(y,y_ref);
    end
    
    %Compare Between The Ground Truth and the Projections Tensor (with the noise) without CP.
    y = double(TensorProjNoise(:,:,proj_num));
    y_ref = double(TensorProj(:,:,proj_num));
    RMSE_NoCp(proj_num) = rmse(y,y_ref);
    PSNR_NoCp(proj_num) = psnr(y,y_ref);
    [SSIM_NoCp(proj_num), ~] = ssim(y,y_ref);
    
    % create RMSE/PSNR/SSIM VECTOR (without CP) for the graph of each proj
    RMSE_NoCp_vec{proj_num,1} = repmat(RMSE_NoCp(proj_num),1,number_of_ranks);
    PSNR_NoCp_vec{proj_num,1} = repmat(PSNR_NoCp(proj_num),1,number_of_ranks);
    SSIM_NoCp_vec{proj_num,1} = repmat(SSIM_NoCp(proj_num),1,number_of_ranks);
end

% Compare full image - Ground Truth and all the Projections with several CP RANK.
for rank = 1:50
    RMSE_vec_ALL(rank) = rmse(TensorProjNoiseCP_ALL{1,rank},TensorProj_ALL{1,rank});
    PSNR_vec_ALL(rank) = psnr(TensorProjNoiseCP_ALL{1,rank},TensorProj_ALL{1,rank});
    [SSIM_vec_ALL(rank), ~] = ssim(TensorProjNoiseCP_ALL{1,rank},TensorProj_ALL{1,rank});
end

%Calc the RMSE/PSNR/SSIM of ALL proj
RMSE_ALL_NoCP = rmse(double(TensorProjNoise),double(TensorProj));
PSNR_ALL_NoCP = psnr(double(TensorProjNoise),double(TensorProj));
[SSIM_ALL_NoCP, ~] = ssim(double(TensorProjNoise),double(TensorProj));

% create RMSE/PSNR/SSIM VECTOR (without CP) for the graph of ALL
RMSE_vec_NoCp_ALL = repmat(RMSE_ALL_NoCP,1,50);
PSNR_vec_NoCp_ALL = repmat(PSNR_ALL_NoCP,1,50);
SSIM_vec_NoCp_ALL = repmat(SSIM_ALL_NoCP,1,50);

%% finding the max of the PSNR and the min of the RMSE for each projection
for proj_num = 1:number_of_projections
    [minRMSEval(proj_num),minRMSEindex(proj_num)] = min(RMSE_mat(:,proj_num));
    [maxPSNRval(proj_num),maxPSNRindex(proj_num)] = max(PSNR_mat(:,proj_num));
    [maxSSIMval(proj_num),maxSSIMindex(proj_num)] = max(SSIM_mat(:,proj_num));
end

% finding the max of the PSNR and the min of the RMSE for ALL projection
[minRMSEval_ALL,minRMSEindex_ALL] = min(RMSE_vec_ALL(:));
[maxPSNRval_ALL,maxPSNRindex_ALL] = max(PSNR_vec_ALL(:));
[maxSSIMval_ALL,maxSSIMindex_ALL] = max(SSIM_vec_ALL(:));

%% Load data
load('P:\Files\ALL_16_8');

%% Create Graph of the Compare Between The Ground Truth and the each Projection with several CP RANK.
RANK = 1:number_of_ranks;
figure(9);
% first projection
subplot(4,1,1);
plot(RANK,RMSE_mat(:,1),'o');
title(['Compare Between The Ground Truth and the ', '\bf FIRST ', '\rm of the Projection with several CP RANK.']);
xlim([0 50]);
grid on;
hold on;
plot(RANK,PSNR_mat(:,1),'o');
hold on;
plot(RANK,SSIM_mat(:,1),'o');
% mark the best RMSE/PSNR/SSIM
plot(minRMSEindex(1) , minRMSEval(1),'--gs');
hold on;
plot(maxPSNRindex(1) , maxPSNRval(1),'--gs'); 
hold on;
plot(maxSSIMindex(1) , maxSSIMval(1),'--gs'); 
hold on;
% adding plots of Reference (without CP)
plot(1:number_of_ranks,RMSE_NoCp_vec{1,1});
hold on;
plot(1:number_of_ranks,PSNR_NoCp_vec{1,1});
hold on;
plot(1:number_of_ranks,SSIM_NoCp_vec{1,1});
hold on;
legend('RMSE','PSNR','SSIM',['RMSE_BEST= (',num2str(minRMSEindex(1)),' , ',num2str(minRMSEval(1)),')'],['PSNR_BEST= (',num2str(maxPSNRindex(1)),' , ',num2str(maxPSNRval(1)),')'],['SSIM_BEST= (',num2str(maxSSIMindex(1)),' , ',num2str(maxSSIMval(1)),')'],...
    ['RMSE_NoCp = ',num2str(RMSE_NoCp(1))],['PSNR_NoCp = ',num2str(PSNR_NoCp(1))],['SSIM_NoCp = ',num2str(SSIM_NoCp(1))]);

% second projection
subplot(4,1,2);
plot(RANK,RMSE_mat(:,2),'o');
title(['Compare Between The Ground Truth and the ', '\bf SECOND ', '\rm of the Projection with several CP RANK.']);
xlim([0 50]);
grid on;
hold on;
plot(RANK,PSNR_mat(:,2),'o');
hold on;
plot(RANK,SSIM_mat(:,2),'o');
% mark the best RMSE/PSNR
plot(minRMSEindex(2) , minRMSEval(2),'--gs');
plot(maxPSNRindex(2) , maxPSNRval(2),'--gs');
plot(maxSSIMindex(2) , maxSSIMval(2),'--gs'); 
% adding plots of Reference (without CP)
plot(1:number_of_ranks,RMSE_NoCp_vec{2,1});
hold on;
plot(1:number_of_ranks,PSNR_NoCp_vec{2,1});
hold on;
plot(1:number_of_ranks,SSIM_NoCp_vec{2,1});
hold on;
legend('RMSE','PSNR','SSIM',['RMSE_BEST= (',num2str(minRMSEindex(2)),' , ',num2str(minRMSEval(2)),')'],['PSNR_BEST= (',num2str(maxPSNRindex(2)),' , ',num2str(maxPSNRval(2)),')'],['SSIM_BEST= (',num2str(maxSSIMindex(2)),' , ',num2str(maxSSIMval(2)),')'],...
    ['RMSE_NoCp = ',num2str(RMSE_NoCp(2))],['PSNR_NoCp = ',num2str(PSNR_NoCp(2))],['SSIM_NoCp = ',num2str(SSIM_NoCp(2))]);


%third projection
subplot(4,1,3);
plot(RANK,RMSE_mat(:,3),'o');
title(['Compare Between The Ground Truth and the ', '\bf THIRD ', '\rm of the Projection with several CP RANK.']);
xlim([0 50]);
grid on;
hold on;
plot(RANK,PSNR_mat(:,3),'o');
hold on;
plot(RANK,SSIM_mat(:,3),'o');
% mark the best RMSE/PSNR
plot(minRMSEindex(3) , minRMSEval(3),'--gs');
plot(maxPSNRindex(3) , maxPSNRval(3),'--gs'); 
plot(maxSSIMindex(3) , maxSSIMval(3),'--gs'); 
% adding plots of Reference (without CP)
plot(1:number_of_ranks,RMSE_NoCp_vec{3,1});
hold on;
plot(1:number_of_ranks,PSNR_NoCp_vec{3,1});
hold on;
plot(1:number_of_ranks,SSIM_NoCp_vec{3,1});
hold on;
legend('RMSE','PSNR','SSIM',['RMSE_BEST= (',num2str(minRMSEindex(3)),' , ',num2str(minRMSEval(3)),')'],['PSNR_BEST= (',num2str(maxPSNRindex(3)),' , ',num2str(maxPSNRval(3)),')'],['SSIM_BEST= (',num2str(maxSSIMindex(3)),' , ',num2str(maxSSIMval(3)),')'],...
    ['RMSE_NoCp = ',num2str(RMSE_NoCp(3))],['PSNR_NoCp = ',num2str(PSNR_NoCp(3))],['SSIM_NoCp = ',num2str(SSIM_NoCp(3))]);

%ALL projection
subplot(4,1,4);
plot(1:50,RMSE_vec_ALL,'o');
title(['Compare Between The Ground Truth and ', '\bf ALL ', '\rm of the Projection with several CP RANK.']);
xlim([0 50]);
grid on;
hold on;
plot(1:50,PSNR_vec_ALL,'o');
hold on;
plot(1:50,SSIM_vec_ALL,'o');
% mark the best RMSE/PSNR/SSIM  
plot(minRMSEindex_ALL , minRMSEval_ALL,'--gs');
hold on;
plot(maxPSNRindex_ALL , maxPSNRval_ALL,'--gs'); 
hold on;
plot(maxSSIMindex_ALL , maxSSIMval_ALL,'--gs'); 
hold on;
%adding plots of Reference(ALL) (without CP)
plot(1:50,RMSE_vec_NoCp_ALL);
hold on;
plot(1:50,PSNR_vec_NoCp_ALL);
hold on;
plot(1:50,SSIM_vec_NoCp_ALL);
hold on;
legend('RMSE','PSNR','SSIM',['RMSE_BEST= (',num2str(minRMSEindex_ALL),' , ',num2str(minRMSEval_ALL),')'],...
    ['PSNR_BEST= (',num2str(maxPSNRindex_ALL),' , ',num2str(maxPSNRval_ALL),')'],...
    ['SSIM_BEST= (',num2str(maxSSIMindex_ALL),' , ',num2str(maxSSIMval_ALL),')'],...
    ['RMSE_NoCp = ',num2str(RMSE_ALL_NoCP)],['PSNR_NoCp = ',num2str(PSNR_ALL_NoCP)],['SSIM_NoCp = ',num2str(SSIM_ALL_NoCP)]);

% Label: ALL 16_8.mat

