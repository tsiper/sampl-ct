%% Exploring the properties of limited angle systems
Initialize;
tic;
N = 32;
K = 32;
M = 32;

Sigma = 1e-3;

ThetaMax = 120;

%% Loading the phantom and generating projections
x0 = LoadPhantom(N);
theta = (0:1/M:1-1/M)*ThetaMax;
Al = paralleltomo(N,theta,K,N);

y0 = reshape(Al*x0(:),[K,M]);
yn = GaussianNoise(y0,Sigma);

%% Exploring some properties of Al
RankAl = rank(full(Al));

%% Attempting naive approach for solving
x1 = vec2im(Al \ yn(:));

%% Plotting
figure;
subplot(221);imagesc(x0);
subplot(222);imagesc(y0);
subplot(223);imagesc(x1);
subplot(224);imagesc(yn);
toc;