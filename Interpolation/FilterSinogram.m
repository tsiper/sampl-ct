function [ y_filt] = FilterSinogram( y)
%FILTERSINOGRAM Receives the sinogram y and filters it assuming a range of 180
%degrees

[m,n] = size(y);
Nfft = 2^(nextpow2(max(m,2*n)+2));

% Nfft = max(m,2*n);

y_pad = [y,flipud(y)];
y_fft = fftshift(fft2(y_pad,Nfft,Nfft));

% [M,N] = size(y_fft);

H_bowtie = BowTie(Nfft,Nfft,0.1);
% Guassian filtering
h_filt = fspecial('gaussian',floor([Nfft/20 Nfft/20]),Nfft/60);
H_bowtie = imfilter(H_bowtie,h_filt);

% H_bowtie = ones(size(y_fft));
% H_bowtie = ones(size(y_rd_fft));
y_fft_filt = y_fft.*(H_bowtie);

y_filt = real(ifft2(ifftshift(y_fft_filt)));
y_filt = y_filt(1:m,1:2*n);

y_filt = [flipud(y_filt(:,end/2+1:ceil(3*end/4))),y_filt(:,ceil(end/4)+1:end/2)];

% Projecting to only positive values
% y_filt(y_filt < 0) = 0;

end

