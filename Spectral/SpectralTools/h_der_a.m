function [ H_der ] = h_der_a( a, Spectra, MassAttenCoeff  )

I = length(a);
data_size = size(a{1});
K = length(Spectra);
for n=1:I
    a_vec(:,n) =a{n}(:); % size: thicknesses X materials
end

denominator = cell(K,1);
nominator = cell(I,K);
H_der = zeros(I,K,data_size(1), data_size(2));

% calculates h derivative for each spectrum and each material:
for k = 1:K
    f_vec = [];
    for n =1:I
        f_vec(n,:) = MassAttenCoeff{n,k}; % size: materials X energies (of spectrum No. k)
    end
    transmission = exp(-a_vec*f_vec);% size: thicknesses X energies
    denominator{k} = reshape(sum(transmission*diag(Spectra{k}.Intensity),2),data_size)+eps;
    for n=1:I
        nominator{n,k} = reshape(sum(transmission*diag(f_vec(n,:).*Spectra{k}.Intensity),2),data_size);
        H_der(n,k,:,:) = (nominator{n,k}./denominator{k});
    end
end


