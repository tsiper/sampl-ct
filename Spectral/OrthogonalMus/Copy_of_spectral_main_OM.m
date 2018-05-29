%% Spectral CT
% Project by Gil Ben Ari and Rotem Zamir
clear; clc;
% to run script you need to add all the CT folder to the path.

%% Loading simulation parameters 
Copy_of_spectral_load_params_OM;

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

%% Calculating orthogonal mus

% MuSpectrum = XrayTubeSpectrumTasmip(max(SpectralParams.Vtube),'number_spectrum','clipzeros');
MuSpectrum = XrayTubeSpectrumTasmip(OM_Vtube,'number_spectrum','clipzeros');
% interpolating for more energies:
x_temp = MuSpectrum.energies;
MuSpectrum.energies = linspace(x_temp(1),x_temp(end),SpectralParams.InterpolationRes);
Mu = zeros(SpectralParams.InterpolationRes,I);

for m=1:I
    Mu(:,m) = getMassAttenCoeff(ChosenMaterials{m}, MuSpectrum.energies);
end
[ort_mus,invP_operator, P_operator, eigs,U] = OrgthogonalMus(Mu');
if ~ReconstructionParams.toOM
    invP_operator = eye(I);
    P_operator = eye(I);
end

%%
figure;
subplot(211);
plot(MuSpectrum.energies,Mu(:,1),'b');
hold on
plot(MuSpectrum.energies,Mu(:,2),'y');
plot(MuSpectrum.energies,Mu(:,3),'g');
xlabel('energies[keV]');
subplot(212);
plot(MuSpectrum.energies,ort_mus(:,1),'k');
hold on
plot(MuSpectrum.energies,ort_mus(:,2),'m');
plot(MuSpectrum.energies,ort_mus(:,3),'r');
xlabel('energies[keV]');

%% Spectral Reconstruction - using FISTA
%Calculate Mass attenuation coefficient for chosen materials
MassAttenCoeff = cell(I,SpectralParams.K ); %Mass attenuation coefficient

for k =1:SpectralParams.K
    temp_mu = zeros(length(Spectra{k}.energies),I);
    for n=1:I
        temp_mu(:,n) = getMassAttenCoeff(ChosenMaterials{n}, Spectra{k}.energies);
    end
    ort_mu = temp_mu*invP_operator;
    for n=1:I
        MassAttenCoeff{n,k} = ort_mu(:,n);
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
    p_input{n} = GaussianNoise(projections{n},ReconstructionParams.noiseLevel);
%     p_input{n} = PoissonNoise(projections{n},35);
end

column_stack_w = zeros(ReconstructionParams.PhantomRes^2,I);
for m=1:I
    column_stack_w(:,m) = phantom{m}(:);
end
column_stack_ort_w = (P_operator*column_stack_w')';

% run recontruction
w_true = cell(I,1);
for m=1:I
    w_true{m} = reshape(column_stack_ort_w(:,m),...
        [ReconstructionParams.PhantomRes,ReconstructionParams.PhantomRes]);
end
[a_hat, f_vals, psnr_vals] = spectral_FISTA_reconstruction_OM(p_input,h,h_grad,a_input,...
    ReconstructionParams,w_true,invP_operator,P_operator);
%% Save
% save(['.\Spectral\results\',save_name,'.mat'])




