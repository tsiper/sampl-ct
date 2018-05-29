function [ ppSinogram ] = InterpolateWindow( Sinogram, N, theta, WindowSize, B, R, W )
% Interpolates the sinogram using the ideal windowed kernel, when given the

global DebugFlag;

% If no B,R,W are given we assume defaults
if nargin < 4
    Mp = 2*N+2;
    R = 1/pi/sqrt(2);
    W = pi*size(Sinogram,1);
    B = floor(Mp/8);
    WindowSize = 8;
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

%% Getting the original samples

% The measurements

A = @(x) real(ifft2(fft2(x).*abs(fft2(H_window))))/numel(x);

%% Filtering with L2 norm
lambda = 1e-2;

% Solving for the optimal solution in the frequency domain
Q = abs(fft2(H_window))/numel(H_window);
b = (fft2(Sinogram));
b_hat = (Q.*b) ./ (Q.^2+lambda);
% b_hat = 
FiltSinogram = real(ifft2(b_hat));

%% Filtering with TV
% lambda = 5e-2;
% L = 50;
% iters = 50;
% TV_Iters = 30;
% % Normalizing the input to the range [0,1]
% Amplitude = range(Sinogram(:));
% Sinogram_Norm = Sinogram/Amplitude;
% FiltSinogram = FISTA_TV(Sinogram_Norm,A,A,zeros(size(Sinogram)),lambda,L,iters,TV_Iters);
% FiltSinogram = FiltSinogram * Amplitude;

%% Commencing Calculation
% Preallocation
ppSinogram  = zeros(2*N+1 , 2*N+2 );

% Starting the waitbar
cpb = ConsoleProgressBar();
% Starting main iteration of MFISTA
cpb.start(); cpb.setText('Subspace Interpolation Algorithm');
% Running the interpolation steps
for l=1:size(ppSinogram,2)
    for k=1:size(ppSinogram,1)

        if (theta_pp(l) >= -45) && (theta_pp(l) < 45)
            K = k-1-(size(ppSinogram,1))/2;
            T = abs(cosd(theta_pp(l)));
        else
            K = k-1-(size(ppSinogram,1)+2*(theta_pp(l)-45)/45)/2;
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
        for j=  (-delta:delta) + min_theta_idx
            for i=(-delta:delta) + min_t_idx
                if (i>0)&&(i<=length(d_t))&&(j>0)&&(j<=length(d_theta))
                    % Distance from center point
                    Rad = sqrt(d_t(i)^2+(d_theta(j)/pi)^2);
                    % Implementing Hamming window
                    Window = 0;
                    if Rad<(WindowSize/2/N)
                        Window = 0.54 + 0.46*cos(2*pi*Rad/WindowSize*N);
                    end
                    % The interpolation kernel 
                    % TODO Precalculate the Sinogram Kernel points
                    h = T*Window*SinogramKernel(d_theta(j),d_t(i),B,R,W)/numel(FiltSinogram);
                    % Performing the interpolation
                    ppSinogram(k,l) = ppSinogram(k,l) + h*FiltSinogram(i,j);
                end
            end
        end
    end
    
    cpb.setValue(l/size(ppSinogram,2));
    
end
cpb.stop();

end
