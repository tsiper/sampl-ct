function [ y_denoise ] = DenoiseSinogram( y,B,R,W,K )
%DENOISESINOGRAM This function receives a sinogram y, and denoises it using
%frequency convolution with the bowtie filter
% 
% % B         = floor(Mp/8);
% B         = floor(Mp/20);
% % R         = 1/pi/2/sqrt(2);
% 
% R        = 1/pi/4;  
% W         = (pi)*(2*N+1);

% Default values for B,R and W

global DebugFlag;

[m,n] = size(y);

% Making the image square by padding with zeros
if m<n
    y_pad = padcols(y,n);
elseif m>n
    y_pad = padcols(y',m)';
end

if nargin==1
    W = pi * size(y_pad,1);
    R = 1/pi/2;
    B = size(y_pad,2)/3;
    K = 20;
end



% 
theta = linspace(-pi/2,pi/2,max(m,n));
t     = linspace(-1/2,1/2,max(m,n));
h     = SinogramKernel( theta,t,B,R,2*W )/n^2/2;

% Using a big window to avoid some problems
w     = HammingWindow(t,theta,K*2*pi/max(m,n));
% w     = ones(length(t),length(theta));

% H = BowTie(max(m,n),max(m,n),0.1);
% h = real(ifft2(ifftshift(H)));

if DebugFlag, ShowImage(h.*w); end

y_denoise = FreqConv2(y_pad,h.*w);

% Getting back to normal size
if m<n
    y_denoise = trimcols(y_denoise,m);
elseif m>n
    y_denoise = trimcols(y_denoise',n)';
end

% Normalizing DC value
y_denoise = y_denoise/(sum(vec(h.*w)));

end

