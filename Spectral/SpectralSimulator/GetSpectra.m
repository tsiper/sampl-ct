function [ Spectra ] = GetSpectra( ScanType, Params, ChosenMaterials, base )
% All Spectra are calculated based on the assumption of TASMIP

%scantype defines the type of the scan for the simulation. It can be:
%'switching' - used mainly for DECT, where the tube's voltage changes (usually 2 different
%voltages). This option simulates the needed spectra.
%'Thresholding' - used for spectral scan, where the detectors are energy
%sensitive. This option simulates the thresholded bins' spectra.

% Params is a struct containing the scan's parameters, depends on the scan type:
% K - number of spectra.
% Vtube - the voltage/s of the tube, in Kev. When in switching, it is a
% vector with K elements.
% N0 - Number of photons emmited (dosage)
% Thresholds - needed for 'Thresholding' scan type. A vector containing the
% threshold energies for calculating the bins.
% InterpolationRes - number of energies samples in the spectra.
% extSpectra - a cell array of the external spectra. containing the fields:
% energies, specnum and Nphotons.

Spectra = cell(Params.K,1);
BodyMaterials = {'adipose', 'blood', 'bone_compact', 'bone_cortical', 'brain', 'lung',...
    'muscle_skeletal', 'muscle_striated', 'skin', 'soft_tissue', 'water'};

switch ScanType
    case 'Thresholding'
        ThresholdSpectra = XrayTubeSpectrumTasmip(Params.Vtube,'number_spectrum','clipzeros');
        % interpolating for more energies:
        x_temp = ThresholdSpectra.energies;
        y_temp = ThresholdSpectra.specnum;
        ThresholdSpectra.energies = linspace(x_temp(1),x_temp(end),Params.InterpolationRes);
        ThresholdSpectra.specnum = interp1(x_temp,y_temp,ThresholdSpectra.energies);
        ThresholdSpectra.specnum_N0 = round(ThresholdSpectra.specnum/ThresholdSpectra.Nphotons*Params.N0);
        
        if Params.ToOptimizeThresholds           
            for n=1:length(ChosenMaterials)
                 ThresholdSpectra.mus(:,n)=getMassAttenCoeff(ChosenMaterials{n},ThresholdSpectra.energies);
            end  
            a_coor = base/2;
            Transmission_E = exp(-ThresholdSpectra.mus*a_coor(:));
            attenuated_spectrum = ThresholdSpectra.specnum_N0(:).*Transmission_E(:);
            specnum = attenuated_spectrum;
            
            if sum(specnum) < 1000 % make sure sum(specnum) big enough so can do integer calcs below
                specnum = 1000*specnum/sum(specnum);
            end
            integral = cumsum(specnum);
            idx = zeros(Params.K-1,1);
            for k = 1:(Params.K-1)
                idx(k) = find(integral< floor(k*integral(end)/Params.K),1,'last');
            end
            Params.Threshold = [ThresholdSpectra.energies(idx), Params.Vtube];
        end
        
        LowThreshIDX = 1;
        for n =1:Params.K
            energydiff = abs(ThresholdSpectra.energies-Params.Threshold(n));
            HighThreshIDX = find(energydiff==min(energydiff));
            Spectra{n}.Intensity = ThresholdSpectra.specnum_N0(LowThreshIDX:HighThreshIDX);
            Spectra{n}.energies = ThresholdSpectra.energies(LowThreshIDX:HighThreshIDX);
            LowThreshIDX = HighThreshIDX+1;
        end
    case 'switching'
        for n =1:Params.K
            SwitchSpectra = XrayTubeSpectrumTasmip(Params.Vtube(n),'number_spectrum','clipzeros');
            % interpolating for more energies:
            x_temp = SwitchSpectra.energies;
            y_temp = SwitchSpectra.specnum;
            SwitchSpectra.energies = linspace(x_temp(1),x_temp(end),Params.InterpolationRes);
            SwitchSpectra.specnum = interp1(x_temp,y_temp,SwitchSpectra.energies);
            Spectra{n}.Intensity = round(SwitchSpectra.specnum/SwitchSpectra.Nphotons*Params.N0);
            Spectra{n}.energies = SwitchSpectra.energies;
        end
    case 'ExtSpectra'
        for n =1:Params.K
            % interpolating for more energies:
            spectrum = Params.extSpectra{n};
            x_temp = spectrum.energies;
            y_temp = spectrum.specnum;
            extSpectra.energies = linspace(x_temp(1),x_temp(end),Params.InterpolationRes);
            extSpectra.specnum = interp1(x_temp,y_temp,extSpectra.energies);
            Spectra{n}.Intensity = round(extSpectra.specnum/extSpectra.Nphotons*Params.N0);
            Spectra{n}.energies = extSpectra.energies;
        end
end

end

