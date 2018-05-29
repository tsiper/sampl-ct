function [ I_hat ] = F2( I )
%F2 Receives an image I and makes a 2D padded dft on it using fft

I_hat = F1( F1( I.' ).' );

end