% Spectral CT
% Project by Gil Ben Ari and Rotem Zamir
Initialize;
% to run script you need to add all the CT folder to the path.

%% Loading simulation parameters
alvarez_params;

%% Loading phantom
load AlvarezPhantom.mat 
RGBphantom  =  alvarezPhantom; %decimate

phantom = cell(I,1);
for ii=1:I
    phantom{ii} = imresize( ReconstructionParams.base(ii)*RGBphantom(:,:,ii),[ReconstructionParams.PhantomRes ReconstructionParams.PhantomRes]);
end

%% Creating spectra
% calling GetSpectra or loading external spectra data
% Spectra = GetSpectra(scan_type, SpectralParams, ChosenMaterials, ReconstructionParams.base);

%% Creating the scan of the object (simulating a spectral scan)
% [ Spectra, projections, RadonMaterials, ChosenMaterials ] = scan( Spectra, phantom, ChosenMaterials );

%% Spectral Reconstruction - using FISTA
%Calculate Mass attenuation coefficient for chosen materials
% MassAttenCoeff = cell(I,SpectralParams.K ); %Mass attenuation coefficient
% for n=1:I
%     for k =1:SpectralParams.K
%         MassAttenCoeff{n,k} = getMassAttenCoeff(ChosenMaterials{n}, Spectra{k}.energies);
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spectrum_data = XrayTubeSpectrumTasmip(SpectralParams.Vtube,'number_spectrum','clipzeros');
% interpolating for more energies:
x_temp = spectrum_data.energies;
y_temp = spectrum_data.specnum;
spectrum_data.energies = linspace(x_temp(1),x_temp(end),SpectralParams.InterpolationRes);
spectrum_data.specnum = interp1(x_temp,y_temp,spectrum_data.energies);
spectrum_data.specnum_N0 = round(spectrum_data.specnum/spectrum_data.Nphotons*SpectralParams.N0);
%% Calculating coefficients of chosen materials at specified energies:
for n=1:I+1
        spectrum_data.mus(:, n) = getMassAttenCoeff(ChosenMaterials{n}, spectrum_data.energies);
end
spectrum_data.mus(:,3) = (1-iodine_frac)*spectrum_data.mus(:,3)+iodine_frac*spectrum_data.mus(:,4);
spectrum_data.mus(:,4) = [];
Spectra = cell(SpectralParams.K,1);
if SpectralParams.ToOptimizeThresholds
    a_coor = calibration_region/2;
    Transmission_E = exp(-spectrum_data.mus*a_coor(:));
    attenuated_spectrum = spectrum_data.specnum_N0(:).*Transmission_E(:);
    specnum = attenuated_spectrum;
    
    if sum(specnum) < 1000 % make sure sum(specnum) big enough so can do integer calcs below
        specnum = 1000*specnum/sum(specnum);
    end
    integral = cumsum(specnum);
    idx = zeros(SpectralParams.K-1,1);
    for k = 1:(SpectralParams.K-1)
        idx(k) = find(integral< floor(k*integral(end)/SpectralParams.K),1,'last');
    end
    SpectralParams.Threshold = [spectrum_data.energies(idx), SpectralParams.Vtube];
end
LowThreshIDX = 1;
for n =1:SpectralParams.K
    energydiff = abs(spectrum_data.energies-SpectralParams.Threshold(n));
    HighThreshIDX = find(energydiff==min(energydiff));
    Spectra{n}.Intensity = spectrum_data.specnum_N0(LowThreshIDX:HighThreshIDX);
    Spectra{n}.energies = spectrum_data.energies(LowThreshIDX:HighThreshIDX);
    LowThreshIDX = HighThreshIDX+1;
end
%% Creating the scan of the object (simulating a spectral scan)
% [ Spectra, projections, RadonMaterials, ChosenMaterials ] = scan( Spectra, phantom, ChosenMaterials );
RadonMaterials = cell(I,1);
for n=1:I
    RadonMaterials{n} = Hpp(phantom{n});
end
K = SpectralParams.K;
MassAttenCoeff = cell(I+1,K);
for k =1:K
    for n=1:I+1  
        MassAttenCoeff{n,k} = getMassAttenCoeff(ChosenMaterials{n}, Spectra{k}.energies);
    end
    MassAttenCoeff{3,k} = (1-iodine_frac)*MassAttenCoeff{3,k}+iodine_frac*MassAttenCoeff{4,k};
    MassAttenCoeff{4,k} = [];  
end


projections = h_a(RadonMaterials, Spectra, MassAttenCoeff );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% run recontruction
w_true = phantom;
[a_hat, f_vals, psnr_vals] = spectral_FISTA_reconstruction_HtH(p_input,h,h_grad,a_input,...
    ReconstructionParams,w_true);
figure;
%% Save
save([save_name,'.mat'])

%%
% post processing with TV
D = DecOperator(2,'uniform');
[ HtH] = BuildAtA( ReconstructionParams.PhantomRes ,2,'uniform' );
for m = 1:I
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

figure; imagesc(w_k{1});
figure; imagesc(w_k{2});


