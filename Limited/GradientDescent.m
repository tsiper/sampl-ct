function [ ErrPlot, f_k] = GradientDescent( N, K, M, ThetaRange, Sigma, f_k, f0 )

% N = The resolution of the image 
% K = number of detectors.
% M =  The number angles/projections
% ThetaRange
% Sigma
Theta = (0:1/M:1-1/M)*ThetaRange + (180-ThetaRange)/2;


A = paralleltomo(N, Theta, K);

Sinogram = A*f0(:);
Sinogram_Noise = GaussianNoise(Sinogram,Sigma);

mu = 2e-5;
f_k = f_k - mu*A'*(A*f_k-Sinogram_Noise);
ErrPlot = norm(A*f_k-Sinogram,'fro');

    
end

