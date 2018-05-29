function [ ppSinogram ] = InterpolateFast( Sinogram, N, theta, WindowSize, B, R, W , lambda )
% Interpolates the sinogram using the ideal windowed kernel, when given the

global DebugFlag;

% Default values if not all are given
if nargin < 4
    Mp = 2*N+2;
    R = 1/pi/sqrt(2);
    W = pi*size(Sinogram,1);
    B = floor(Mp/8);
    WindowSize = 8;
    lambda = 1e-2;
end

% The Pseudo-Polar scanning angles
l = -N/2:N/2;
theta_pp = [atan(2*l/N),fliplr(pi/2-atan(2*l/N))]*180/pi;

% Building the ideal Filtering kernel
t_vec     = linspace(-1/2,1/2,size(Sinogram,1));%-1/(2*size(y_rd_pad,1));
theta_vec = linspace(-1/2,1/2,size(Sinogram,2))*(max(theta)-min(theta));
H_full    = SinogramKernel(deg2rad(theta_vec),t_vec,B,R,W);
% Window    = HammingWindow(t_vec,theta_vec/(max(theta)-min(theta)),2*WindowSize/size(Sinogram,1));
Window    = HammingWindow(t_vec,theta_vec/180,2*WindowSize/size(Sinogram,1));
H_window  = H_full.*Window;

% Normalizing the center value to be 1
H_window = H_window./(sum(H_window(:)))*numel(H_window);

if DebugFlag
    ShowImage(H_window);
end

%% Solving for the optimal solution in the frequency domain + L2 Norm
Q = abs(fft2(H_window))/numel(H_window);
b = (fft2(Sinogram));
b_hat = (Q.*b) ./ (Q.^2+lambda);
FiltSinogram = real(ifft2(b_hat));

%% Commencing Calculation
% Preallocation
ppSinogram  = zeros(2*N+1 , 2*N+2 );
[Mpp,Npp] = size(ppSinogram);

% Running the interpolation steps
parfor l=1:Npp
    for k=1:Mpp

        if (theta_pp(l) >= -45) && (theta_pp(l) < 45)
            K = k-1-(Mpp)/2;
            T = abs(cosd(theta_pp(l)));
        else
            K = k-1-(Mpp+2*(theta_pp(l)-45)/45)/2;
            T = abs(sind(theta_pp(l)));
        end

        I =  (1:size(FiltSinogram,1)) - 1 - (size(FiltSinogram,1)+1)/2;
        d_theta = deg2rad( theta_pp(l)-theta );
        d_t     = (K*T-I)/size(FiltSinogram,1);
        
        % Another small correction found experimentally
        d_t = d_t - abs( mod(theta_pp(l)+45,90) - 45 )/180/size(FiltSinogram,1);
        
        %% Fast window design
        % Finding the correct range of points
        [~,min_theta_idx] = min(abs(d_theta));
        [~,min_t_idx    ] = min(abs(d_t));
        delta = ceil(WindowSize);

        % Iterating over all the relevant points with contribution
        i =  (-delta:delta) + min_t_idx;
        i((i<1)|(i>length(d_t))) = [];
        j =  (-delta:delta) + min_theta_idx;
        j((j<1)|(j>length(d_theta))) = [];
        
        % Building the window function
        [d_Theta,d_T] = meshgrid(d_theta(j),d_t(i));
        Rad = sqrt((d_T).^2+(d_Theta/pi).^2);
        Window = 0.54 + 0.46*cos(2*pi*Rad/WindowSize*N);
        Window(Rad>=(WindowSize/2/N)) = 0;
        
        % Computing the relevant kernel
        h = T*Window.*SinogramKernel(d_theta(j),d_t(i),B,R,W)/numel(FiltSinogram);

        % Performing the interpolation for a specific point
        ppSinogram(k,l) = vec(h)'*vec(FiltSinogram(i,j));

    end
       
end

end
