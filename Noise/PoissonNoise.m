function [ Yn ] = PoissonNoise( Y,SNR )
%POISSONNOISE Receives an array Y, and adds poisson noise such that the SNR is
%equal to SNR we assume that Y is positive

% Normalizing and forcing positivity
Y(Y<0) = 0;
MaxY = max(Y(:));
Ynorm = Y/MaxY;

% Calculating the noise power according to SNR
Py = norm(Ynorm(:)).^2;
Pn = Py*10^(-SNR/10);

% Finding the Poisson distribution with a variance of Pn
% lambda = 1e3*Ynorm./mean(Ynorm(:))/(Pn);
% R  = poissrnd(lambda);
% Renormalizing
% Yn = R./max(R(:))*MaxY;

R = sqrt(Ynorm).*randn(size(Ynorm));

% Finding the constant for correct SNR
A = 10^(-SNR/10)*Py/norm(R(:))^2;

Yn = Ynorm + sqrt(A)*R;
Yn(Yn<0) = 0;
% Yn = abs(Yn);
Yn = Yn*MaxY;
% Yn = 1e8*imnoise(1e-11*Y,'poisson');

end

