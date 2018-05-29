%% Spectral CT
% Project by Gil Ben Ari and Rotem Zamir

clear
close all;
clc;
% to run script you need to add all the CT folder to the path.

%% Spectra Initialization
% Number of spectra involved in each scan
SpectralParams.K = 3;
% Amount of photons emitted, proportional to dosage
SpectralParams.N0 = 10^20;
% The tube voltage, if we are on "switching" it should have length(K) elements
% else (thresholding) just the center frequency of the single emitted spectrum
SpectralParams.Vtube = 140;
% The bin right limits markers, the first bin starts at 0
SpectralParams.Threshold = [60 ,100, 140];
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
% % FIXME: A normalizing vector, trying to equalize energies between materials
base = [0.125 1 20]./255;
% ChosenMaterials = { 'blood','adipose', 'water'}; %choose the materials you wish for the basis
% % FIXME: A normalizing vector, trying to equalize energies between materials
% base = [1,1,1];
% ChosenMaterials = { 'brain','bone_compact'}; %choose the materials you wish for the basis
% % FIXME: A normalizing vector, trying to equalize energies between materials
% base = [1.052,1.145];
% ChosenMaterials = {'adipose', 'water'}; %choose the materials you wish for the basis
% % FIXME: A normalizing vector, trying to equalize energies between materials
% base = [1,1];

%% Loading phantom
% chosen phantom can be either: 'circles', 'more_circles', 'paint',
% 'spectral_shepp', 'spectral_zubal'.
chosen_phantom = 'circles';
RGBphantom = LoadPhantom(PhantomRes,chosen_phantom);
% RGBphantom = zeros(PhantomRes,PhantomRes,I);
% RGBphantom(:,:,1) = rgb2gray(im2double(imresize(imread('brain_map.png'),[PhantomRes,PhantomRes])));
% RGBphantom(:,:,2) = rgb2gray(im2double(imresize(imread('bone_map.png'),[PhantomRes,PhantomRes])));
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
% simulate the corresponding attenuation (the spectrum is interpulated
% through GetSpectra()).

% calling GetSpectra or loading external spectra data
Spectra = GetSpectra('Thresholding', SpectralParams);

%%

%NEED TO DO: Normalize phantom levels to actual range of thicknesses (graylevel = 255);

%% Creating the scan of the object (simulating a spectral scan)
[ Spectra, projections, RadonMaterials, ChosenMaterials ] = scan( Spectra, phantom, ChosenMaterials );

%% Spectral Reconstruction - using FISTA
global DebugFlag PlotFlag;
DebugFlag=1;
PlotFlag=1;

%Parameters for reconstruction
lambda = 0.01; %lambda for TV step
L =800*SpectralParams.K; %Lipschitz coeeficient. optimal 90 with calibration region , 800 for distinct materials,0.5adipose&water
iters = 10000; %Maximal number of iterations of the minimization algorithm
TV_iters = 20; %Number of iterations in TV
L_coeffs = ones(1,SpectralParams.K); %P metric coefficients, used for scaling
toTV =0; %Flag stating whether to perform TV or not
toStop = 0;
MassAttenCoeff = cell(I,SpectralParams.K ); %Mass attenuation coefficient

%Calculate Mass attenuation coefficient for chosen materials
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
    p_input{n} = projections{n};
%     p_input{n} = GaussianNoise(projections{n},0.005);
%     p_input{n} = PoissonNoise(projections{n},35);
end
% run recontruction
w_true = phantom;
[a_hat, f_vals, psnr_vals] = MFISTA_recontruction(p_input,h,h_grad,a_input,...
    lambda,L,iters,TV_iters,RadonMaterials,PhantomRes,L_coeffs,toTV,toStop,base,w_true);
figure;
%% Save
save('init_rnd0dot02_else02.mat')
%%
MyPP = @(x) real(invF1(App(x))); %Pseudo Polar transformation
MyAdjPP = @(y) real(App_T(M(F1(y))));
AppAdjApp = @(x) MyAdjPP(MyPP(x));
MyInversePP = @(y) conjgrad( AppAdjApp , MyAdjPP(y), zeros(PhantomRes,PhantomRes) , 6);

alpha_hat = zeros(PhantomRes,PhantomRes,I);
alpha_hat(:,:,1) = MyInversePP(a_hat{1})./base(1);
alpha_hat(:,:,2) = MyInversePP(a_hat{2})./base(2);
alpha_hat(:,:,3) = MyInversePP(a_hat{3})./base(3);

true_alpha = zeros(PhantomRes,PhantomRes,I);
true_alpha(:,:,1) = MyInversePP(RadonMaterials{1})./base(1);
true_alpha(:,:,2) = MyInversePP(RadonMaterials{2})./base(2);
true_alpha(:,:,3) = MyInversePP(RadonMaterials{3})./base(3);
%% TV reconstruction
lambda_rec = 0.01;
TV_Iters_rec = 20;
alpha_hat_TV = zeros(PhantomRes,PhantomRes,I);
for m = 1:I
    
    % Finding the norm of the spatial image for material m
    range_k = range(vec(alpha_hat(:,:,m)));
    min_k = min(vec(alpha_hat(:,:,m)));
    arg_k_norm = (alpha_hat(:,:,m)-min_k)/range_k;
    
    %perform FGP
    w_k = FGP(arg_k_norm, 2*lambda_rec, TV_Iters_rec);
    
    % Normalizing back the image to its original intensity
    alpha_hat_TV(:,:,m) = w_k*range_k+min_k;
end
%% plotting
figure()
subplot(231);imagesc(a_hat{1}); colorbar('southoutside');
subplot(232);imagesc(a_hat{2}); colorbar('southoutside');
subplot(233);imagesc(a_hat{3}); colorbar('southoutside');

subplot(234);imshow(alpha_hat_TV,[]);colorbar('southoutside');title('reconstructed phantom');
colorbar('southoutside');
subplot(235);imshow(true_alpha,[]);colorbar('southoutside');title('original phantom');
diff = abs(true_alpha-alpha_hat_TV);
subplot(236);imshow(diff,[]);colorbar('southoutside');title('difference');

% subplot(234);imshow(alpha_hat,[]);colorbar('southoutside');title('reconstructed phantom');
% colorbar('southoutside');
% subplot(235);imshow(true_alpha,[]);colorbar('southoutside');title('original phantom');
% diff = abs(true_alpha-alpha_hat);
% subplot(236);imshow(diff,[]);colorbar('southoutside');title('difference');











