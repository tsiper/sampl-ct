%% Spectral CT
% Project by Gil Ben Ari and Rotem Zamir
clear; clc;
% to run script you need to add all the CT folder to the path.

%% Loading simulation parameters 
spectral_load_params_6;

%% Loading phantom
RGBphantom = LoadPhantom(ReconstructionParams.PhantomRes,chosen_phantom);
phantom = cell(IReal,1);
for ii=1:IReal
    phantom{ii} = ReconstructionParams.baseReal(ii)*RGBphantom(:,:,ii);
end

%% Creating spectra
% calling GetSpectra or loading external spectra data
Spectra = GetSpectra(scan_type, SpectralParams, ChosenMaterialsReal, ReconstructionParams.baseReal);

%% Creating the scan of the object (simulating a spectral scan)
[ Spectra, projections, RadonMaterials, ChosenMaterialsReal ] = scan( Spectra, phantom, ChosenMaterialsReal );

%% Spectral Reconstruction - using FISTA
%Calculate Mass attenuation coefficient for chosen materials
MassAttenCoeff = cell(IRecon,SpectralParams.K ); %Mass attenuation coefficient
for n=1:IRecon
    for k =1:SpectralParams.K
        MassAttenCoeff{n,k} = getMassAttenCoeff(ChosenMaterialsRecon{n}, Spectra{k}.energies);
    end
end

%Defining projection functions
h = @(a) (h_a(a, Spectra, MassAttenCoeff));
h_grad = @(a) h_der_a(a, Spectra, MassAttenCoeff);

% MFISTA input
a_input = cell(IRecon,1);
p_input = cell(SpectralParams.K,1);
for m = 1:IRecon
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

% run recontruction
w_true = phantom([2,3]);
[a_hat, f_vals, psnr_vals] = spectral_FISTA_reconstruction_HtH(p_input,h,h_grad,a_input,...
                                ReconstructionParams,w_true);
figure;
%% Save
save(['.\Spectral\results\',save_name,'.mat'])

%%
% post processing with TV
D = DecOperator(2,'uniform');
[ HtH] = BuildAtA( ReconstructionParams.PhantomRes ,2,'uniform' );
   for m = 1:IRecon
        b = Hpp_T(D(a_hat{m}));
        temp_w = conjgrad(HtH,b,zeros(size(b)),100,1e-40);
        temp_w(temp_w<0)=0;
        
        min_w = min(temp_w(:));
        range_w = range(temp_w(:));
        norm_w_k = (temp_w-min_w)/range_w;
        norm_x_k = FGP( norm_w_k ,0.07,20);
        w_k{m} = norm_x_k*range_w+min_w;
        z_k{m}  = Hpp(w_k{m});
    end

figure; imagesc(w_true{1});
figure; imagesc(w_true{2});

