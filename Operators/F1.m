function [ g ] = F1( f )
%F1 Performs 1D dft for the signal f and outputs the results g

% Length of columns
N = size(f,1);
% Working on a row vector too
if N==1,  N = size(f,2); end

% Computing the FFT
g = fftshift(fft(f),1);

% Fixing the modulation problem
% If N is dividable by 4 or (N-1) is dividable by 4
if  (~mod(N,4))||(~mod(N-1,4))
    g(2:2:end,:) = -g(2:2:end,:);
else
    g(1:2:end,:) = -g(1:2:end,:);
end