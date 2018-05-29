%% Spectral CT
% Based on Alvarez's work. Project by Gil Ben Ari and Rotem Zamir

clear
close all;
clc;
% addpath(genpath('P:\'))

%% Needed - Based on materials, Choosing A space boundaries

%% Constants initialization

K = 4; %number of energy bins
N0 = 10^20; %number of total photons emmited at source
I = 3; %number of materials
Vtube = 120; %Tube voltage, in kev
theta = 0:2:179; %angles for scanning
T = 100; %number of detectors at each angle
Y = 64; %phantom resolution
interpolation_res = 1000; %energy resolution
gray_level = 255; %gray levels of the phantom\
calibration_res = 10; %calibration space resolution
calibration_region =  [0.125 1 20];%boundries of materials in calibration
plot_flag = 0;
save_name = 'Result291216.mat';

% make sure number of angles and Y hold 
%% choosing materials and phantom
% materials for MusFromMaterials func (will be inside BodyMaterials):
% adipose, blood, bone_compact, bone_cortical, brain, lung, muscle_skeletal,
% muscle_striated, skin, soft_tissue, water

% materials for XrayMu func will be string with Chemical formula or simple
% atomic notation

BodyMaterials = {'adipose', 'blood', 'bone_compact', 'bone_cortical', 'brain', 'lung',...
    'muscle_skeletal', 'muscle_striated', 'skin', 'soft_tissue', 'water'};
ChosenMaterials = { 'I','Ca', 'soft_tissue'}; %choose the materials you wish for the basis
% chosen phantom can be either: 'circles', 'more_circles', 'paint',
% 'spectral_shepp', 'spectral_zubal'.
chosen_phantom = 'circles';
%% Extracting Spectral Data
% Choosing the desired materials and creating their spectral information
    % Calculating TASMIP spectrum:
spectrum_data = XrayTubeSpectrumTasmip(Vtube,'number_spectrum','clipzeros');
    % interpolating for more energies:
x_temp = spectrum_data.energies;
y_temp = spectrum_data.specnum;
spectrum_data.energies = linspace(x_temp(1),x_temp(end),interpolation_res);
spectrum_data.specnum = interp1(x_temp,y_temp,spectrum_data.energies);
spectrum_data.specnum_N0 = round(spectrum_data.specnum/spectrum_data.Nphotons*N0);
%% Calculating coefficients of chosen materials at specified energies:
for n=1:I
    if any(cellfun(@(x) isequal(x, ChosenMaterials{n}), BodyMaterials))
        spectrum_data.mus(:, n) = MusFromMaterials(ChosenMaterials(n), spectrum_data.energies);
    else
        spectrum_data.mus(:, n) = XrayMu(ChosenMaterials{n},spectrum_data.energies);
    end
end
%% calibrator data 

a = zeros(I, calibration_res);
calibration_data.a_calib_region = calibration_region;   %boundaries of calibration region
Material_norm = gray_level./calibration_data.a_calib_region;

for n=1:I
    a(n,:) = linspace(0,calibration_data.a_calib_region(n),calibration_res);
end
calibration_data.As = nMeshGrid(transpose(a));
%% energy bin thresholds calculation
[spectrum_data.idx_threshold, atten_specnum] = OptimNbinsThresholds(spectrum_data,K, 'a_space', calibration_data.a_calib_region/2);
[calibration_data.p_k, calibration_data.ns, calibration_data.N0] = PHAbinCounts(calibration_data.As,spectrum_data,spectrum_data.idx_threshold);
[Q_MLE, a_Error, calibration_data.AsInitFit] = CalibrationAlgorithm(calibration_data);
%% Phantom scanning (for us, craeting Atable from phantom)
phantom = LoadPhantom(Y,chosen_phantom);
MyRadon = @(X) radon(X,theta,T);
scan_data.As = zeros(T,length(theta),I);
for ii =1:I
    scan_data.As(:,:,ii)=MyRadon(phantom(:,:,ii)/Material_norm(ii));
end

if plot_flag
    figure();
    imshow(phantom,[])
    title('Original phantom')
end

%% object scan
[scan_data.p_k, scan_data.ns, scan_data.N0] = PHAbinCounts(scan_data.As,spectrum_data,spectrum_data.idx_threshold);

if plot_flag
    figure();
    p_k_rec = zeros(Y,Y,K);
    for k=1:K
        p_k_rec(:,:,k) = iradon(scan_data.p_k(:,:,k),theta,'linear','Ram-Lak',1,Y);
        subplot(K,1,k);
        imshow(p_k_rec(:,:,k),[]);
        colorbar;
    end
end

%% saving scan_data

save(['P:\git\CT\spectral\algorithm',save_name],'calibration_data','scan_data','spectrum_data','Q_MLE','Material_norm');

%% reconstruction
p_k_size = size(scan_data.p_k);
p_k_vector = reshape(scan_data.p_k, p_k_size(1)*p_k_size(2),p_k_size(3));
A_hat_0_vec = p_k_vector*Q_MLE';

% Do interpolation of the a_error by griddata

A_hat=A_hat_0_vec;
alpha_hat = zeros(Y,Y,I);
for ii=1:I
    alpha_hat(:,:,ii) = iradon(Material_norm(ii)*reshape(A_hat(:,ii),[T , size(theta,2)]),theta,'linear','Ram-Lak',1,Y);
end
 %% Result plot
if plot_flag %needs to be adjusted to new dimensions
    figure();
    quiver3(calibration_data.As(:,1),calibration_data.As(:,2),calibration_data.As(:,3),a_Error(:,1),a_Error(:,2),a_Error(:,3));
    scan_error = A_hat_0-scan_data.As;
    figure();
    quiver3(scan_data.As(:,1),scan_data.As(:,2),scan_data.As(:,3),scan_error(:,1),scan_error(:,2),scan_error(:,3));

    figure();
    scatter3(scan_data.As(:,1),scan_data.As(:,2),scan_data.As(:,3),'b.');
    hold on
    scatter3(A_hat(:,1), A_hat(:,2), A_hat(:,3))
    hold off

    figure();
    scatter3(calibration_data.As(:,1),calibration_data.As(:,2),calibration_data.As(:,3),'b.');
    hold on
    scatter3(calibration_data.AsInitFit(:,1), calibration_data.AsInitFit(:,2), calibration_data.AsInitFit(:,3))
    hold off
end

h1 = figure();
subplot(3,2,1)
imshow((alpha_hat(:,:,1)) , []);
title(['Reconstruction of the color red, material ', ChosenMaterials{1}, ' calib region ' num2str(calibration_data.a_calib_region(1))]);
colorbar;
subplot(3,2,2)
imshow((phantom(:,:,1)) , []);
title('Original color of red');
colorbar;
subplot(3,2,3)
imshow((alpha_hat(:,:,2)) , []);
title(['Reconstruction of the color green, material ', ChosenMaterials{2}, ' calib region ' num2str(calibration_data.a_calib_region(2))]); 
colorbar;
subplot(3,2,4)
imshow((phantom(:,:,2)) , []);
title('Original color of green');
colorbar;
subplot(3,2,5)
imshow((alpha_hat(:,:,3)) , []);
title(['Reconstruction of the color blue, material ', ChosenMaterials{3}, ' calib region ' num2str(calibration_data.a_calib_region(3))]); 
colorbar;
subplot(3,2,6)
imshow((phantom(:,:,3)) , []);
title('Original color of blue');
colorbar;

figure();imshow((alpha_hat) , []);
title('Reconstruction of the body image');
colorbar;

%%
if plot_flag
    %plot tube spectrum
    figure()
    plot(spectrum_data.energies,spectrum_data.specnum_N0);
    xlabel('Energies [Kev]');
    ylabel('Intensity [a.u.]');
    
    %plot energy thresholds
    figure()
    hold on
    plot(spectrum_data.energies,atten_specnum);
    y_lims = get(gca, 'ylim');
    Energy1 = spectrum_data.energies(spectrum_data.idx_threshold(1));
    Energy2 = spectrum_data.energies(spectrum_data.idx_threshold(2));
    Energy3 = spectrum_data.energies(spectrum_data.idx_threshold(3));
    line([Energy1 Energy1], y_lims, 'color', 'r', 'LineStyle', ':');
    line([Energy2 Energy2], y_lims, 'color', 'r', 'LineStyle', ':');
    line([Energy3 Energy3], y_lims, 'color', 'r', 'LineStyle', ':');
    hold off
    xlabel('Energies [Kev]');
    ylabel('Intensity')
    
    %plot attenuation coefficients
    figure()
    hold on
    for n=1:I
        plot(spectrum_data.energies,spectrum_data.mus(:,n));
    end 
end