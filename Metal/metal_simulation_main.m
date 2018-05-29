%% Spectral CT
% Project by Gil Ben Ari and Rotem Zamir
Initialize;
% to run script you need to add all the CT folder to the path.

%% Loading simulation parameters 
load_params_metal_simulation;

%% Loading phantom

RGBphantom = LoadPhantom(ReconstructionParams.PhantomRes,chosen_phantom);
phantom = cell(I,1);
for ii=1:I
    phantom{ii} = ReconstructionParams.base(ii)*RGBphantom(:,:,ii);
end
[ReconstructionParams, I, ChosenMaterials, phantom] = AddMetal(ReconstructionParams, I, ChosenMaterials, phantom);
%% Creating spectra
% calling GetSpectra or loading external spectra data
Spectra = GetSpectra(scan_type, SpectralParams, ChosenMaterials, ReconstructionParams.base);

%% Creating the scan of the object (simulating a spectral scan)

dtheta = 0.3;
theta = 0:dtheta:(180-dtheta);
[ Spectra, projections, RadonMaterials, ChosenMaterials ] = scan_metal( Spectra, phantom, ChosenMaterials,theta );

%% Reconstruction using FBP

reconstructed_phantom = iradon(projections{1},theta);
contaminated_sinogram = radon(reconstructed_phantom,theta);

figure;
subplot(2,2,1);
imagesc(projections{1});
colorbar;
colormap(gray);
title('Original Sinogram');

subplot(2,2,2);
imagesc(reconstructed_phantom);
colorbar;
colormap(gray);
title('FBP');

subplot(2,2,3);
imagesc(contaminated_sinogram);
colorbar;
colormap(gray);
title('Contaminated Sinogram - radon(FBP)');

subplot(2,2,4);
imagesc(contaminated_sinogram(2:end-1,:)-projections{1});
colorbar;
colormap(gray);
title('Artifact Sinogram');


