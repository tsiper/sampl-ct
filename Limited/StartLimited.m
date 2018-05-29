% This is a limited test file created by Dima, Rotem and Shahar

%% Initializing
Initialize();

%% Setting up parameters

N = 256; % The resolution of the image
M = 3; % The number of projections
Sigma = 0.01; % The noise variance

% The angular range
ThetaRange = 120;
Theta = (0:1/M:1-1/M)*ThetaRange + (180-ThetaRange)/2;

%% Loading the phantom
x0 = LoadPhantom(N,'brain');
imshow(x0)
%% Defining the system matrix and projecting the phantom
[A, b, x] = paralleltomo(N,Theta,N,N-1);
A_op  = @(x) reshape(A*vec(x),[N,M]);
At_op = @(y) vec2im(A'*vec(y));

y0  = A_op(x0);
y  = GaussianNoise(y0,Sigma);

%% Traditional reconstruction
x_bp  = At_op(y);
x_fbp = vec2im(fbp(A,vec(y),Theta));

%% Advanced reconstruction
x_start = zeros(size(x0));
lambda = 0.1;
iters = 500; TV_Iters = 30;
L = 1000*N;
x_tv  = FISTA_TV( y,A_op,At_op,x_start,lambda,L,iters,TV_Iters,x0);
%% Real inverse
x_inv = zeros(size(x0));
if N<=64
    x_inv = vec2im(A\vec(y));
end

%% Plotting
figure;
subplot(231);imagesc(x0);   axis square;colormap bone;
subplot(232);imagesc(y);    axis square;colormap bone;
subplot(233);imagesc(x_bp); axis square;colormap bone;
subplot(234);imagesc(x_fbp);axis square;colormap bone;
subplot(235);imagesc(x_tv); axis square;colormap bone;
subplot(236);imagesc(x_inv);axis square;colormap bone;

%% Comparison
CompareImages(x_tv,x0);   