%% Spectral CT
% Project by Gil Ben Ari and Rotem Zamir
clc;
clear all;
Initialize;
% to run script you need to add all the CT folder to the path.

%% Loading simulation parameters 
spectral_load_params_2;

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
[ Spectra, projections, RadonMaterials, ChosenMaterials ] = scan( Spectra, phantom, ChosenMaterials );

%% Spectral Reconstruction - using FISTA
%Calculate Mass attenuation coefficient for chosen materials
MassAttenCoeff = cell(I,SpectralParams.K ); %Mass attenuation coefficient
for n=1:I
    for k =1:SpectralParams.K
        MassAttenCoeff{n,k} = getMassAttenCoeff(ChosenMaterials{n}, Spectra{k}.energies);
    end
end

%Defining projection functions
h = @(a) (h_a(a, Spectra, MassAttenCoeff));
h_grad = @(a) h_der_a(a, Spectra, MassAttenCoeff);

% MFISTA input
a_input = cell(I,1);
p_input = cell(SpectralParams.K,1);
for m = 1:I
    %choose initial guess for the algorithm:
%     a_input{m} = 0.0*RadonMaterials{m}+0.02*randn(size(RadonMaterials{m}));
    a_input{m}=zeros(size(RadonMaterials{m}));
% a_input{m} = RadonMaterials{m};
end

for n=1:SpectralParams.K
    % choose the projections input:
%     p_input{n} = projections{n};
%     p_input{n} = GaussianNoise(projections{n},ReconstructionParams.noiseLevel);
%     p_input{n} = PoissonNoise(projections{n},35);
      p_input{n} = projections{n} + ReconstructionParams.noiseLevel.*range(projections{n}).*randn(size(projections{n})); %simple noise

end

%%  Clean The Noise with Tensor CP Decompositions
[ p_input_CP , p_input_noised , CP_Rank ] = Denoising_CP( projections , p_input , SpectralParams.K ,...
                                                                            'AlreadyNoised' ,0,-1,220); 

%% CP RECONSTRUCTION
% run recontruction
w_true = phantom;
[a_hat_CP, f_vals_CP, psnr_vals_CP,w_best_CP] = spectral_FISTA_reconstruction_HtH_Tensor(p_input_CP,h,h_grad,a_input,...
                                ReconstructionParams,w_true);
%%
figure(15), hold on; suptitle('Reconstruction with and without CP decomposition. simple noise. sigma = 0.05 ' );
subplot(2,3,1);
imagesc(RGBphantom);
title('The Ground True');

subplot(2,3,6);
loglog(1:length(psnr_vals_CP) , psnr_vals_CP);
title(['The Best Psnr (with CP): ', num2str(max(psnr_vals_CP))]);
xlim([0 length(psnr_vals_CP)]);
grid on;

alpha_hat(:,:,1) = w_best_CP{1}/ReconstructionParams.base(1);
alpha_hat(:,:,2) = w_best_CP{2}/ReconstructionParams.base(2);
alpha_hat(:,:,3) = w_best_CP{3}/ReconstructionParams.base(3);
alpha_hat = (alpha_hat-min(alpha_hat(:)))/range(alpha_hat(:));
subplot(2,3,3);
imagesc(alpha_hat);
title('The Reconstructed image with CP');
hold on;

%second figure
figure(30), hold on; suptitle('Reconstruction with and without CP decomposition. simple noise. sigma = 0.05 ' );
subplot(3,1,2);
plot(1:length(psnr_vals_CP) , psnr_vals_CP);
title(['The Best Psnr (with CP): ', num2str(max(psnr_vals_CP))]);
xlim([0 length(psnr_vals_CP)]);
grid on;

%third figure
figure(45), hold on; suptitle('Reconstruction with and without CP decomposition. simple noise. sigma = 0.05 ' );
subplot(1,3,2);
loglog(1:length(psnr_vals_CP) , psnr_vals_CP);
title(['The Best Psnr (with CP): ', num2str(max(psnr_vals_CP))]);
xlim([0 length(psnr_vals_CP)]);
grid on;

%% NO CP RECONSTRUCTION
% run recontruction
w_true = phantom;
[a_hat_NOCP, f_vals_NOCP, psnr_vals_NOCP , w_best_NOCP] = spectral_FISTA_reconstruction_HtH_Tensor(p_input,h,h_grad,a_input,...
                                ReconstructionParams,w_true);
%%
figure(15), hold on;

subplot(2,3,5);
loglog(1:length(psnr_vals_NOCP) , psnr_vals_NOCP);
title(['The Best Psnr (without CP): ', num2str(max(psnr_vals_NOCP))]);
xlim([0 length(psnr_vals_NOCP)]);
grid on;

alpha_hat(:,:,1) = w_best_NOCP{1}/ReconstructionParams.base(1);
alpha_hat(:,:,2) = w_best_NOCP{2}/ReconstructionParams.base(2);
alpha_hat(:,:,3) = w_best_NOCP{3}/ReconstructionParams.base(3);
alpha_hat = (alpha_hat-min(alpha_hat(:)))/range(alpha_hat(:));
subplot(2,3,2);
imagesc(alpha_hat);
title('The Reconstructed image without CP');
hold off;

%second figure
figure(30), hold on;
subplot(3,1,3);
plot(1:length(psnr_vals_NOCP) , psnr_vals_NOCP);
title(['The Best Psnr (without CP): ', num2str(max(psnr_vals_NOCP))]);
xlim([0 length(psnr_vals_NOCP)]);
grid on;

%third figure
figure(45), hold on;
subplot(1,3,3);
loglog(1:length(psnr_vals_NOCP) , psnr_vals_NOCP);
title(['The Best Psnr (without CP): ', num2str(max(psnr_vals_NOCP))]);
xlim([0 length(psnr_vals_NOCP)]);
grid on;
%% Save

%post processing with TV
D = DecOperator(2,'uniform');
[ HtH] = BuildAtA( ReconstructionParams.PhantomRes ,2,'uniform' );
   for m = 1:I
        b = Hpp_T(D(a_hat_NOCP{m}));
        temp_w = conjgrad(HtH,b,zeros(size(b)),100,1e-40);
        temp_w(temp_w<0)=0;
        
        min_w = min(temp_w(:));
        range_w = range(temp_w(:));
        norm_w_k = (temp_w-min_w)/range_w;
        norm_x_k = FGP( norm_w_k ,0.07,20);
        w_k{m} = norm_x_k*range_w+min_w;
        z_k{m}  = Hpp(w_k{m});
    end

figure; imagesc(w_k{1});
figure; imagesc(w_k{2});

