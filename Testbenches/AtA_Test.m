%% Intitializing the workspace
Initialize();
DebugFlag = 0;

%% Run Parameters
N = 256;
DecFactor = 3; % We have 2*N/DecFactor projection angles
Sigma = 0.02;
LoopAll = 0;   % Calculate all the different matrices, overriding DecFactor

%% Loading the phantom
x   = LoadPhantom(N,'zubal');

% Designing the decimation factor, to determine how many viewing angles we have
D = DecOperator(DecFactor,'uniform');

% Bulding all the combos and overwriting them accordingly
if LoopAll
    for N = [64,128,256]
        for DecFactor = [2,3,4,5,6,8]
            % Getting the AtA operator
            [ AtA ] = BuildAtA( N ,DecFactor,'uniform' , 1 );
        end
    end
else
    [ AtA ] = BuildAtA( N ,DecFactor,'uniform' );
end
%% Scanning using Radon

dt = 1;                             % How many degrees we jump in each scan
theta = -90:dt:180;                 % The angular range - this should be bigger than -45:135
y   = radon(x,theta);               % The radon projection
[y_n,SNR] = GaussianNoise(y,Sigma); % Adding noise

%% Preparing the data for Pseudo-Polar solver

% Resampling the noisy measurements to Pseudo-Polar grid
y_pp = InterpolateSinogram(y_n,2*N,theta,'spline');

% If we want to generate directly the PP project than:
x_pad = padarray(x,[N/2,N/2]);
y_pp2 = Rpp(x_pad);

% Generating the measurement vector with all the operators such that
%    y=A*x   ==>   b = At*y = AtA x
b = real(App_T(D(M(F1(y_pp)))));
% Trimming unnecessary zeros
b = b(N/2+1:3*N/2,N/2+1:3*N/2);

% Finding the optimal Lipshitz constant of our AtA system
L = PowerMethod(AtA,N);

%% Running Fista for solving
tic;
Iters = 100; % 20 should be enough - in some cases 8-12 is also enough

x_tv = FISTA_TV_PP(b,AtA,zeros(size(x)),L,Iters);

toc;
double(x_tv);

% Analysing the results
CompareImages(x_tv,x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OLD CODE for debug purposes - Please disregard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x1 = padarray(x1,[N/2,N/2]);
% %%
% % x1 = real(conv2(x,h,'same'));
% % x_pad = repmat(x,[3 3]);
% 
% % %%
% % [Ipp1,Ipp2] = OptimizedPPFT(x);
% % x3 = real(OptimizedadjPPFT(D[(M(Ipp1)),M(Ipp2))));

%% Calculating speed gain
iters = 5;
tic;
for i=1:iters
    x1 = AtA(x);
end
t_conv = toc;
tic;
for i=1:iters
    x1 = App_T(D(M(App(x))));
end
t_old = toc;

time_improvement = t_old / t_conv;



%%
% tic;
% x2 = real(App_T(M(D(App(padarray(x,[N/2 N/2]))))));
% x2 = x2(N/2+1:3*N/2,N/2+1:3*N/2);
% 
% % x2 = real(App_T(M(D(App(x)))));
% 
% toc;
% % figure;
% % imagesc(real(x2-x3));

% %%
% figure;
% subplot(231);imagesc(x);  title('Ground Truth');
% subplot(233);imagesc(x1); title('Frequency Convolution');
% subplot(234);imagesc(x2); title('$App^{T} M D App$');
% subplot(232);imagesc(x2-x1); title('Diff'); colorbar;
% subplot(235);imagesc(x_tv);