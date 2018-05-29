function     [OutputPhantom1NFFT,achivedSNR] = AddNoiseToData(OriginalOutputPhantom1NFFT,DesiredSNR);
NoiseSeed = 1;

signal = OriginalOutputPhantom1NFFT(:,3) + 1i.*OriginalOutputPhantom1NFFT(:,4);
Psignal = mean(abs(signal).^2);
disp(['Signal power is ',num2str(10*log10(Psignal)),' dB']);

% rng(NoiseSeed);
noise = (randn(size(signal)) + 1i.*randn(size(signal))).*sqrt(10^(-DesiredSNR/10)*Psignal/2);
Pnoise = mean(abs(noise).^2);
disp(['Noise power is ',num2str(db(Pnoise)/2),' dB']);
achivedSNR = 10*log10(Psignal/Pnoise);

disp(['Achived SNR is ',num2str(achivedSNR),' dB']);
output = signal + noise;

OutputPhantom1NFFT(:,1) = OriginalOutputPhantom1NFFT(:,1);
OutputPhantom1NFFT(:,2) = OriginalOutputPhantom1NFFT(:,2);
OutputPhantom1NFFT(:,3) = real(output);
OutputPhantom1NFFT(:,4) = imag(output);

end

