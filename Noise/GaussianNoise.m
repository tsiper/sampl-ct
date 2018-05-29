function [ x_noise,SNR ] = GaussianNoise( x,sigma )
%GUASSIANNOISE Adding strictly positive Gaussian noise, with a std of sigma

% Adding the Gaussian noise
x_noise = x + sigma*range(x(:))*randn(size(x));

% Moving negative values to the positive side - This actually makes the SNR
% worse
x_noise(x_noise<0) = -x_noise(x_noise<0);

SNR = -20*log10(norm(x_noise-x,'fro')/norm(x,'fro'));


end

