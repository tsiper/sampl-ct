function [ Spectra, projections, a, materials ] = scan( Spectra, phantom, materials )
% Spectra -  will contain the spectra of the spectral scan, without any attenuation.
% phantom  - is a cell containing the density images for each material. Number
% of elements in the cell is the number of materials. All the density
% images should have the same size.
% materials - is a cell, containing the names of the materials

I = length(materials);
K = length(Spectra);
MassAttenCoeff = cell(I,K);
a = cell(I,1);
phantomRes = size(phantom{1},1);

%% Calculating coefficients of chosen materials at spectra energies:
for n=1:I
    for k =1:K
        MassAttenCoeff{n,k} = getMassAttenCoeff(materials{n}, Spectra{k}.energies);
    end
end

%% Scanning using Radon
dt = 1;                             % How many degrees we jump in each scan
theta = -90:dt:180;                 % The angular range - this should be bigger than -45:135

for n=1:I
%     y   = radon(phantom{n},theta);% The radon projection
    % [y_n,SNR] = GaussianNoise(y,Sigma); % Adding noise

    % Resampling the measurements to Pseudo-Polar grid
%     a{n} = InterpolateSinogram(y,2*phantomRes,theta,'spline');
    a{n} = Hpp(phantom{n});
end

%% object scan using h function as described in paper
[projections ] = h_a( a, Spectra,MassAttenCoeff );


