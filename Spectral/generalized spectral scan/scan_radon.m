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

%% Calculating coefficients of chosen materials at spectra energies:
for n=1:I
    for k =1:K
        MassAttenCoeff{n,k} = getMassAttenCoeff(materials{n}, Spectra{k}.energies);
    end
end

%% Phantom scanning in pseudo-polar
MyPP = @(x) real(invF1(App(x))); %Pseudo Polar transformation
for n=1:I
    a{n} = radon(phantom{n});
end 

%% object scan using h function as described in paper
[projections ] = h_a( a, Spectra,MassAttenCoeff );


