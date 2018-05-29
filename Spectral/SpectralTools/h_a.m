function [ projections ] = h_a( a, Spectra,MassAttenCoeff )
% this function calculate the spectal projections p for each spectrum in Spectra, 
% for a set of thicknesses (in a-space). 

I = length(a);
data_size = size(a{1});
K = length(Spectra);
for n=1:I
    a_vec(:,n) =a{n}(:); % size: thicknesses X materials
end

N_detector = cell(K,1);
N_source= cell(K,1);
projections = cell(K,1);

% calculates P for each spectrum:
for k = 1:K
    f_vec = [];
    for n =1:I
        f_vec(n,:) = MassAttenCoeff{n,k}; % size: materials X energies (of spectrum No. k)
    end
    transmission = exp(-a_vec*f_vec);% size: thicknesses X energies
    N_detector{k} = floor(reshape(sum(transmission*diag(Spectra{k}.Intensity),2),data_size));
    N_source{k} = sum(Spectra{k}.Intensity);
    projections{k} = -log(N_detector{k}/N_source{k}+eps);
end
