function [ Denoising_Projections_CP , Denoising_Projections , CP_Rank ] = Denoising_CP( projections , Noised_Projections ,   K , noise_type , sigma , CP_Rank, number_of_ranks) 
 
% projections :  the projection to denoise 
% K : Number of spectra involved in each scan 
% noise_type : The type of noise added to the signal 
% sigma : sigma of noise 
% CP_Rank : the rank for the cp_als decomposition 
%           if if the CP_Rank set to -1 , find the best rank in the range: 
%           1 - "number_of_ranks" 


%% Defining the Tensors from the scan cell data
ProjectionTensor = zeros([size(projections{1}),length(projections)]);
for i=1:length(projections)
    ProjectionTensor(:,:,i) = projections{i};
end

% Defining the Projection Tensor With the Tensor-ToolBox
TensorProj=tensor(ProjectionTensor); %% define the tensor with the toolbox


%% Defining the Noise Tensors

% if the input was noised
ProjectionTensorWithNoise = zeros([size(Noised_Projections{1}),length(Noised_Projections)]);
for i=1:length(projections)
    ProjectionTensorWithNoise(:,:,i) = Noised_Projections{i};
end

%if we want to add noise to the input
if noise_type ~= 'AlreadyNoised'
    % Defining the Projection Tensor With Gussian / Simple Noise 
    ProjectionTensorWithNoise = zeros([size(projections{1}),length(projections)]);
    for i=1:K
        % GaussianNoise
        switch noise_type 
            % GaussianNoise 
            case 'GaussianNoise'         
                [ProjectionTensorWithNoise(:,:,i),~] = GaussianNoise( ProjectionTensor(:,:,i),sigma );

            % simple noise
            case 'SimpleNoise'
             ProjectionTensorWithNoise(:,:,i) = ProjectionTensor(:,:,i) + sigma.*range(ProjectionTensor(:,:,i)).*randn(size(ProjectionTensor(:,:,i)));
        end
    end
end
%define Tensor
TensorProjNoise=tensor(ProjectionTensorWithNoise);

%% The CP Section. if the CP_Rank = -1 , find out the best rank (according to PSNR) and make the CP with that rank.
% if the CP_Rank != -1 , make the CP with CP_Rank

if CP_Rank == -1
    % Collect Data for Compare Between The Ground Truth and the Projections with several CP RANK.
    % Extremly slow script 
    TensorProjNoiseCP_ALL = cell(1,number_of_ranks);
    TensorProj_ALL = repmat({double(TensorProj)},1,number_of_ranks);
    
    % Starting the waitbar
    counter = ConsoleProgressBar();
    % Starting main iteration of MFISTA
    counter.start(); counter.setText('Finding the best CP Decomposition Rank..');
    
    for rank = 1:number_of_ranks
        TensorProjNoiseCP_ALL{1,rank} = double(cp_als(TensorProjNoise,rank));
        % Collect data of the PSNR between the Tensor+Noise+CP and the clean Tensor
        PSNR_vec_ALL(rank) = psnr(TensorProjNoiseCP_ALL{1,rank},TensorProj_ALL{1,rank});
        counter.setValue(rank/number_of_ranks);
    end
    counter.stop();
    
    %find the max PSNR (value and index)
    [maxPSNRval_ALL , CP_Rank] = max(PSNR_vec_ALL(:));
end

%define the tensor after CP with CP_Rank
TensorProjNoiseCP = double(cp_als(TensorProjNoise,CP_Rank));
 
%% Graph Section
%collect data for the comparing
%Calc the /PSNR/ of ALL proj
PSNR_ALL_NoCP = psnr(double(TensorProjNoise),double(TensorProj));
% create /PSNR/ VECTOR (without CP) for the graph of ALL
PSNR_vec_NoCp_ALL = repmat(PSNR_ALL_NoCP,1,number_of_ranks);

% Create Graph of the Compare Between The Ground Truth and the Projections with several CP RANK.
figure(15);

%ALL projection
subplot(2,3,4);
semilogx(1:number_of_ranks,PSNR_vec_ALL);
title('PSNR of the Proj after different decompositions');
xlim([0 number_of_ranks]);
xlabel('Rank');
ylabel('PSNR');
grid on;
hold on;

% mark the best RMSE/PSNR/SSIM  
semilogx(CP_Rank , maxPSNRval_ALL,'--gs'); 
hold on;

%adding plots of Reference(ALL) (without CP)
semilogx(1:number_of_ranks,PSNR_vec_NoCp_ALL);
hold on;

legend('PSNR',['PSNR BEST= (',num2str(CP_Rank),' , ',num2str(maxPSNRval_ALL),')'],...
   ['PSNR NoCp = ',num2str(PSNR_ALL_NoCP)]);

printf(['the CP-rank is ' , num2str(CP_Rank) ])

%% create another graph
% Create Graph of the Compare Between The Ground Truth and the Projections with several CP RANK.
figure(30);

%ALL projection
subplot(3,1,1);
plot(1:number_of_ranks,PSNR_vec_ALL);
title('PSNR of the Proj after different decompositions');
xlim([0 number_of_ranks]);
xlabel('Rank');
ylabel('PSNR');
grid on;
hold on;

% mark the best RMSE/PSNR/SSIM  
plot(CP_Rank , maxPSNRval_ALL,'--gs'); 
hold on;

%adding plots of Reference(ALL) (without CP)
plot(1:number_of_ranks,PSNR_vec_NoCp_ALL);
hold on;

legend('PSNR',['PSNR BEST= (',num2str(CP_Rank),' , ',num2str(maxPSNRval_ALL),')'],...
   ['PSNR NoCp = ',num2str(PSNR_ALL_NoCP)]);

%% create another third graph
% Create Graph of the Compare Between The Ground Truth and the Projections with several CP RANK.
figure(45);

%ALL projection
subplot(1,3,1);
semilogx(1:number_of_ranks,PSNR_vec_ALL);
title('PSNR of the Proj after different decompositions');
xlim([0 number_of_ranks]);
xlabel('Rank');
ylabel('PSNR');
grid on;
hold on;

% mark the best RMSE/PSNR/SSIM  
semilogx(CP_Rank , maxPSNRval_ALL,'--gs'); 
hold on;

%adding plots of Reference(ALL) (without CP)
semilogx(1:number_of_ranks,PSNR_vec_NoCp_ALL);
hold on;

legend('PSNR',['PSNR BEST= (',num2str(CP_Rank),' , ',num2str(maxPSNRval_ALL),')'],...
   ['PSNR NoCp = ',num2str(PSNR_ALL_NoCP)]);



%% Convert the Tensors to cell of matrixes
% Defining the projections cell from the projections Tensor (before CP).
Denoising_Projections_CP =  TensorToMatrixInCell(TensorProjNoiseCP,K);  
Denoising_Projections = TensorToMatrixInCell(TensorProjNoise,K);  
end

%bla