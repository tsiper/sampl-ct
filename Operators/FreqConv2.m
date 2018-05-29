function [ C ] = FreqConv2( A,B )
%FREQCONV Performs frequency domain convolution between A and B


[n1,m1] = size(A);
[n2,m2] = size(B);

if (n1~=m1) || (n2~=m2) || (n1~=n2)
    error('The images A and B must be square and of the same size');
end
% n = n1;

n = ceil(n1/4);

% Padding the arrays
Apad = padarray(A,[n,n]);
Bpad = padarray(B,[n,n]);

% Convlving in frequency domain
Cpad = ifftshift(ifft2(fft2(Apad).*fft2(Bpad)));

C = Cpad(n+1:end-n,n+1:end-n);


end

