function [uniform_k_samples, OutputImage] = FilterInKSpaceandIFFT(c,r_aq,sqrtN,OverGridFactor)

    FOV_idx = ((sqrtN*OverGridFactor)/2-sqrtN/2+1:(sqrtN*OverGridFactor)/2+sqrtN/2);
    d=conv2(r_aq.',r_aq.',c,'same');
    
    IFFT_d = ifftshift(ifft2(fftshift(d)));
    OutputImage = IFFT_d(FOV_idx,FOV_idx);
    OutputImage = real(OutputImage); % Image should be real
    OutputImage = OutputImage.*((sqrtN*OverGridFactor)^2); % Correct gain from oversampling
    
    % Calculate k-space samples on the uniform grid
    OutputImageMasked = zeros(size(IFFT_d));
    OutputImageMasked(FOV_idx,FOV_idx) = OutputImage;
    OutputImageMasked = real(OutputImageMasked);
    uniform_k_samples = ifftshift(fft2(fftshift(OutputImageMasked)));
end

