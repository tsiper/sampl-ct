function [uniform_k_samples, OutputImage] = IFFTandFilterInImageSpace(c,ImgFilter,sqrtN,OverGridFactor)

FOV_idx = ((sqrtN*OverGridFactor)/2-sqrtN/2+1:(sqrtN*OverGridFactor)/2+sqrtN/2);

IFFT_c = ifftshift(ifft2(fftshift(c)));
OutputImage = (sqrtN*OverGridFactor)^2.*IFFT_c.*abs(ImgFilter);

%calculate k-space samples on the uniform grid
uniform_k_samples = ifftshift(fft2(fftshift(real(OutputImage))));

% Crop to FOV and returen real Image
OutputImage = real(OutputImage(FOV_idx,FOV_idx)); 
end

