%% Introduction.
% This is a code-file which allows to find the lowest conditional number
% based on the amount of detectors, angular range and other stuff we will
% think of.

%% Initializing.
Initialize();

%% Playground - A) Wide Theta Range


N = 64; % The resolution of the image + number of detectors.
M = 50; % The number angles/projections

% The angular range
ThetaRange = linspace(10,180, 86); % 18 options.
% Theta = (0:1/M:1-1/M)*ThetaRange + (180-ThetaRange)/2;
Condition_Num = zeros(1,length(ThetaRange));
i=1;

for ThetaRange=linspace(10, 180, 86)
    Theta = (0:1/M:1-1/M)*ThetaRange + (180-ThetaRange)/2;
    [A] = paralleltomo(N,Theta,N,N-1);
    A_full = full(A);
    tic
    Condition_Num(i) = cond(A_full'*A_full);
    toc
    i=i+1;
end

figure(1)
plot(linspace(10,180, 86), Condition_Num)
title('Wide Theta Range')
xlabel('Theta Range')
ylabel('Cond Number')

%% Playground - B) Narrow Theta Range

N = 64; % The resolution of the image + number of detectors.
M = 50; % The number angles/projections

% The angular range
ThetaRange_N = linspace(50,169, 120); % 30 options.
% Theta = (0:1/M:1-1/M)*ThetaRange + (180-ThetaRange)/2;
Condition_Num_N = zeros(1,length(ThetaRange_N));
i=1;

for ThetaRange_N=linspace(50, 169, 120)
    
    Theta = (0:1/M:1-1/M)*ThetaRange_N + (180-ThetaRange_N)/2;
    [A_N] = paralleltomo(N,Theta,N,N-1);
    A_full_N = full(A_N);
    tic
    Condition_Num_N(i) = cond(A_full_N'*A_full_N);
    toc
    i=i+1;
end

figure(2)
plot(linspace(50, 169, 120), Condition_Num_N)
title('Narrow Theta Range')
xlabel('Theta Range')
ylabel('Cond Number')

%% Playground - C) Number of Projections vs. Angular Range

N = 64; % The resolution of the image + number of detectors.
M = linspace(14,64,6); % The number angles/projections
Cond_Num_Mat = zeros(60, 6);

% The angular range
ThetaRange = linspace(50,169, 60); % 30 options.

i=1;
for ThetaRange=linspace(50, 169, 60)
    j=1;
    for M=linspace(14,64,6)
        Theta = (0:1/M:1-1/M)*ThetaRange + (180-ThetaRange)/2;
        [A] = paralleltomo(N,Theta,N,N-1);
        A_full = full(A);
        tic
        Cond_Num_Mat(i,j) = cond(A_full'*A_full);
        toc
        j=j+1;
    end
    i=i+1;
end


figure(3)
[X,Y] = meshgrid(linspace(50, 169, 60),linspace(14,64,6));
surf(X,Y,log(Cond_Num_Mat'))
title('Number of Projections vs. Angular Range')
xlabel('Angular Range')
ylabel('Number of Projections')
zlabel('Logaritmic Condition Number')

% ????? ????? ???? ?? ?? ?????? ???????? ???? ????? ??????? ???? ?????
% ???????, ????? E ????? ????? ?? ???? ?????? ????? ???????, ?????? ?????
% ?? ????? ???, ???? ????? ????? C ???? ????? ?????????. ???? ?????? ?????
% ?? ??? ????? ?????? ?????? ????, ??? ???? ??????? ???, ????? ?????? ?????
% ??? ????. ??? ?? ???? ?? ????? ??? ???? ????? ???????-?????.

%% Playground - D) Verification - Minor Sigma - varying Projections

Sigma = 0.000001; % The noise variance
N = 64; % The resolution of the image 
K = 64; % number of detectors.
ThetaRange = 170;
f0 = LoadPhantom(N);
iters = 100;


% Range- 170;   Projection- 14; 
M_1 = 14; % The number angles/projections
f_k = zeros(size(f0(:)));
figure(4);

ErrPlot = [];
for k = 1:iters
    [ErrPlot(k), f_k] = GradientDescent( N, K, M_1, ThetaRange, Sigma, f_k, f0 ); 
    subplot(321);imagesc(vec2im(f_k));colormap('bone');
    title(['Iter ',num2str(k), '; Minor Sigma; Projection- 14']);
    subplot(322);loglog(ErrPlot);grid on;
    drawnow;
    pause(0.01);
end

% Range- 170;   Projection- 34; 
M_2 = 34; % The number angles/projections
f_k = zeros(size(f0(:)));

ErrPlot = [];
for k = 1:iters
    [ErrPlot(k), f_k] = GradientDescent( N, K, M_2, ThetaRange, Sigma, f_k, f0 ); 
    subplot(323);imagesc(vec2im(f_k));colormap('bone');
    title(['Iter ',num2str(k), '; Minor Sigma; Projection- 34']);
    subplot(324);loglog(ErrPlot);grid on;
    drawnow;
    pause(0.01);
end

% Range- 170;   Projection- 34; 
M_3 = 64; % The number angles/projections
f_k = zeros(size(f0(:)));

ErrPlot = [];
for k = 1:iters
    [ErrPlot(k), f_k] = GradientDescent( N, K, M_3, ThetaRange, Sigma, f_k, f0 ); 
    subplot(325);imagesc(vec2im(f_k));colormap('bone');
    title(['Iter ',num2str(k), '; Minor Sigma; Projection- 64']);
    subplot(326);loglog(ErrPlot);grid on;
    drawnow;
    pause(0.01);
end


%% Playground - E) Verification - Minor Sigma - varying Theta Range

Sigma = 0.000001; % The noise variance
N = 64; % The resolution of the image 
K = 64; % number of detectors.
M = 34;
f0 = LoadPhantom(N);
iters = 100;


% Range- 60;   Projection- 34; 
ThetaRange_1 = 60; % The number angles/projections
f_k = zeros(size(f0(:)));
figure(5);

ErrPlot = [];
for k = 1:iters
    [ErrPlot(k), f_k] = GradientDescent( N, K, M, ThetaRange_1, Sigma, f_k, f0 ); 
    subplot(321);imagesc(vec2im(f_k));colormap('bone');
    title(['Iter ',num2str(k), '; Minor Sigma; ThetaRange- 60']);
    subplot(322);loglog(ErrPlot);grid on;
    drawnow;
    pause(0.01);
end

% Range- 100;   Projection- 34; 
ThetaRange_2 = 100; % The number angles/projections
f_k = zeros(size(f0(:)));

ErrPlot = [];
for k = 1:iters
    [ErrPlot(k), f_k] = GradientDescent( N, K, M, ThetaRange_2, Sigma, f_k, f0 ); 
    subplot(323);imagesc(vec2im(f_k));colormap('bone');
    title(['Iter ',num2str(k), '; Minor Sigma; ThetaRange- 100']);
    subplot(324);loglog(ErrPlot);grid on;
    drawnow;
    pause(0.01);
end

% Range- 160;   Projection- 34; 
ThetaRange_3 = 160; % The number angles/projections
f_k = zeros(size(f0(:)));

ErrPlot = [];
for k = 1:iters
    [ErrPlot(k), f_k] = GradientDescent( N, K, M, ThetaRange_3, Sigma, f_k, f0 ); 
    subplot(325);imagesc(vec2im(f_k));colormap('bone');
    title(['Iter ',num2str(k), '; Minor Sigma; ThetaRange- 160']);
    subplot(326);loglog(ErrPlot);grid on;
    drawnow;
    pause(0.01);
end



%% Playground - F) Number of Projections vs. Number of Detectors

N = 64; % The resolution of the image 
K = linspace(10,60,11); % Number of detectors.
M = linspace(14,64,6); % The number angles/projections
Cond_Num_Mat = zeros(11, 6);
ThetaRange = 170;

i=1;
for K = linspace(10,60,11)
    j=1;
    for M=linspace(14,64,6)
        Theta = (0:1/M:1-1/M)*ThetaRange + (180-ThetaRange)/2;
        [A] = paralleltomo(N,Theta,K,N-1);
        A_full = full(A);
        tic
        Cond_Num_Mat(i,j) = cond(A_full'*A_full);
        toc
        j=j+1;
    end
    i=i+1;
end


figure(6)
[X,Y] = meshgrid(linspace(10,60,11),linspace(14,64,6));
surf(X,Y,log(Cond_Num_Mat'))
title('Number of Projections vs. Number of Detectors')
xlabel('Number of Detectors')
ylabel('Number of Projections')
zlabel('Logaritmic Condition Number')
