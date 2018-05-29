function [ Rsa ] = SinogramSamplingMtx( theta_vec, t_vec,B,R,W,WindowSize )
%BUILDSAMPLINGMTX Summary of this function goes here
%   Detailed explanation goes here

M = length(t_vec);
N = length(theta_vec);


Rsa = zeros(N*M);
% 
% B = length(theta_vec)/10;
% R = 1/pi;
% W = pi*length(t_vec);

% Window = HammingWindow(0:M-1,theta_vec,WindowSize);
Window = hamming(2*WindowSize);
Window = [Window(WindowSize+1:end); zeros(ceil(sqrt(2)*max(N,M)),1)];

for i = 1:N*M
    
    % Setting the k,l indices according to the vectorization
    k = mod((i-1),M)+1;        % Columns of the PP transform
    l = floor((i-1)/M)+1;      % Rows of the PP transform
    
    % Looping over all the columns of App, corresponsing to pixels in the image
    for j=1:N*M
        % Setting the i,j indices according to the vectorization
        m = mod((j-1),M)+1;   % Columns of the original image
        n = floor((j-1)/M)+1; % Rows of the original image
        
        t     = (t_vec(m) - t_vec(k));
        theta = deg2rad(theta_vec(n)-theta_vec(l));
        
        Rsa(i,j) = SinogramKernel(theta,t,B,R,W)*Window(round(sqrt((k-m)^2+(l-n)^2))+1);
%         Rsa(i,j) = SinogramKernel(theta,t,B,R,W);
    end
end

end

