%% Psuedo-Polar reconstruction using iterative steps

%% Initializing and globals
global PlotFlag DebugFlag SaveFlag FIGURE_PATH RandSeed detector_pattern;
FIGURE_PATH = './';

detector_pattern = [];
RandSeed = randi(1e6);
PlotFlag  = 1;
DebugFlag = 1;
SaveFlag  = 0;

%% Initializing
cg_iters = 5;
% size of the image - Always use powers of 2
N = 128;
Mproj  = 2*(N+1);%DecFactor*(N+1)*2;
sigma_x = 0.0;       % Image noise
sigma_y = 0.0;       % Measurement noise
SNR = 300;
PhantomType = 'zubal';

%% Decimation parameters
DecFactor = 4;

%DecMethod = 'uniform';
% DecMethod = 'random';
% DecMethod = 'limited';
DecMethod = 'detector';
%DecMethod = 'pattern algorithm';
%DecMethod = 'detector_mult_theta';
x0 = LoadPhantom(N,PhantomType);
Wx = 1;
Wy = 50;

D = DecOperator(DecFactor,DecMethod);

%% MFISTA Parms
% lambda       = [16 2] * 1e-3;
% lambda       = 24e-3;
lambda       = 1e-8;
L            =    10;
% mfista_iters =   25;
mfista_iters = 100;
TV_Iters     =   10;

%% SFISTA Params
% SFISTA_Params.Lambda       = 50e-3;
% SFISTA_Params.mu0          = 1e-1;
% SFISTA_Params.muf          = 1e-4;
% SFISTA_Params.mu           = SFISTA_Params.mu0;
% % SFISTA_Params.mu           = (1e-3)/SFISTA_Params.Lambda;
% % SFISTA_Params.mu           = 1e-4; 
% SFISTA_Params.gamma        = 10;
% SFISTA_Params.Continuation = 0;
% SFISTA_Params.IterMax      = 200;
% % SFISTA_Params.L            = norm(A, 2)^2 + (Dnorm^2)/SFISTA_Params.mu;        % If D is unitary then norm(D, 2) = 1, I think
% SFISTA_Params.L            = 8;

%% Wavelet operator setting up
% WaveletType =  'Haar', 'Beylkin', 'Coiflet', 'Daubechies',
%                'Symmlet', 'Vaidyanathan','Battle'
% WaveParams.Type     = 'Haar';
% % The size of each filter of the wavelet
% WaveParams.FiltSize = 2;
% % The scale of the wavelet - the number of wavelet levels
% WaveParams.Scale    = log2(N);
% WaveOp = Wavelet(WaveParams.Type,WaveParams.FiltSize,WaveParams.Scale);
% 
% % Weighting -First try is to nullify the first quadrand, related to DC
% WeightMatrix = [0 1; 1 0];
% WeightMatrix = WeightMatrix / norm(WeightMatrix,'fro');
% WaveParams.Weight  = kron( WeightMatrix, ones(2*N) );
% 
% W  = @(t) imresize((WaveOp'*(WaveParams.Weight.*t)),1/4);
% Wt = @(t) WaveParams.Weight.*(WaveOp*imresize(t,4));


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   STARTING THE ALGORITHM
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loading the phantom image
%x0 = LoadPhantom(N,PhantomType);

%x0 = phantom(N);

% Adding noise
x0_noise = x0 + sigma_x * range(x0(:)) * randn(size(x0));

%% FBP and PP calculations
% s       = (-1+4/Mproj):(4/Mproj):1;
theta   = (-1/4:1/Mproj:3/4-1/Mproj)*180;

l = linspace(-N/2,N/2,Mproj/2);
theta_pp = [atan(2*l/N),fliplr(pi/2-atan(2*l/N))]*180/pi;
% theta_pp   = [atan(s),atan(s)+pi/2]*180/pi;
tic;
y_radon = radon(x0_noise,theta);
% y_radon_noise = y_radon + sigma_y*range(y_radon(:))*randn(size(y_radon));
y_radon_noise   = PoissonNoise(y_radon,SNR);
y_pp    = radon(x0_noise,theta_pp);
x_radon = iradon(y_radon_noise,theta);
toc;
x_radon = x_radon(2:end-1,2:end-1);

% X_Radon denoising with TV so we compare apples to apples
x_radon_tv = FGP(x_radon,lambda(end),150);
% y     = A(x0_noise);
 ppSinogram = Radon2pp(y_pp,theta_pp);
% ppSinogram = Radon2ppSinc(y_pp,theta_pp);
% ppSinogram   = InterpolateSpline(y_radon,N,theta,theta_pp);

% This also denoises the sinogram by a moderate ammount
%ppSinogram = Radon2pp2(y_pp,theta_pp);

% ppSinogram = Radon2pp2(y_pp,theta_pp);

% m = size(x0,1)*2+1;
% ppSinogram = (padarray(ppSinogram,[(m-size(ppSinogram,1))/2,0]));
% ppSinogram = FilterSinogram(ppSinogram);
% ppSinogram_Sinc = padarray(ppSinogram_Sinc,[(m-size(ppSinogram_Sinc,1))/2,0]);
ppSinogram = padcols(ppSinogram,2*N+1);
% ppSinogram_noise = ppSinogram + sigma_y*range(ppSinogram(:))*randn(size(ppSinogram));
ppSinogram_noise = PoissonNoise(ppSinogram,SNR);
y = F1(ppSinogram_noise);
% SHAHAR A temporary hack 
% y = A(x0);

y_pp0 = real(invF1(App(x0,Mproj/2)));

PSNR_Sinogram = psnr(ppSinogram_noise,y_pp0);

% y_sinc = F1(ppSinogram_Sinc);
% y_old = F1(ppSinogram_old);
% x0_shift   = [zeros(length(x0),1),x0(:,1:end-1)];


%% Conj Grad Method
cg_iters = 10;
% w = y;

% G = @(x) x;
% Defining the preconditioned operators
AtAop = @(t) App_T((Mnew(App(t,Mproj/2))));
% AtA = @(t) At_op(M(A_op(t)));
AtY = App_T(Mnew(App(x0)));
% AtY_old = App_T(M(y_old));



tic;
x_cg = conjgrad(AtAop,AtY,zeros(size(x0)),cg_iters);
% x_cg_old = conjgrad(AtA,AtY_old,zeros(size(x0)),cg_iters);
toc;

x_cg = vec2im(x_cg);
psnr_cg = psnr(x_cg,x0);
% psnr_cg_old = psnr(x_cg_old,x0_shift)
% ax2(2) = subplot(235);imagesc(abs(y));

%% Solving with MFISTA TV Minimization
% 


% G = @(y) y;

A_Fista  = @(x) (D(Mnew(App(x,Mproj/2))));

if strcmpi(DecMethod,'limited')
    At_Fista = @(y) 1/(1-DecFactor)*App_T(D(y));
else
    At_Fista = @(y) DecFactor*App_T(D(y));
end



% A_Fista  = @(x) D(M(A_op(x)),DecFactor,method);
% At_Fista = @(y) App_T(D(y,DecFactor,method));

% At_Fista = @(y) At_op((D((y),DecFactor,method)));
% y_Fista  = (D(M(App(x0)),DecFactor,method));
y_Fista  = D(Mnew(y));
x_Fista  = zeros(size(x0));
% h = AtA(N,DecFactor,method);
% h = zeros(N);
% x_tv  =  MFISTA(y_Fista,A_Fista,At_Fista,x_Fista,10*lambda,L,mfista_iters/2,TV_Iters,x0);
% x_tv  =  MFISTA(y_Fista,A_Fista,At_Fista,x_Fista,lambda,L,mfista_iters,TV_Iters,x0);

x_tv  =  MFISTA_Bi(y_Fista,A_Fista,At_Fista,x_Fista,lambda,L,mfista_iters,TV_Iters,x0);



% Offset of 1 pixel
% x_tv =  [x_tv(:,2:end), zeros(size(x_tv,1), 1)];
% x_wtv = wMFISTA(y_Fista,A_Fista,At_Fista,zeros(size(x0)),lambda,L,mfista_iters,TV_Iters,Wx,Wy);

% x_tv = x_tv ./ max(x_tv(:));

%% Limited angle frequency boosting
% % x_tv_limited = x_tv;
% % theta_final = 135;
% % theta_start = 135 - 180/DecFactor;
% % H_Wedge     = 1-Wedge(N,N,theta_start,theta_final);
% % 
% % 
% % if strcmpi(method,'limited')
% %     for i=1:MultiLambdaIters
% %         x_tv_fft      = fftshift(fft2(x_tv_limited));
% %         x_tv_fft_filt = x_tv_fft.*(H_Wedge+5)/6;
% %         x_tv_filt     = real(ifft2(ifftshift( x_tv_fft_filt ) ));
% %         x_tv_limited  = MFISTA(y_Fista,A_Fista,At_Fista,x_tv_filt,3e-3,L,5,TV_Iters);
% %     end
% %     
% %     figure;
% %     subplot(221); imagesc(x_tv);
% %     subplot(222); imagesc(x_tv_filt);
% %     subplot(223); imagesc(db(abs(x_tv_fft)));
% %     subplot(224); imagesc(db(abs(x_tv_fft_filt)));
% %     
% % end


%% Iterative Lambdas with warm start
% x_rtv = zeros(size(x0));
% lambda_k = 2*lambda;
% for i=1:MultiLambdaIters
%     if i >  1
%         close;
%     end
% %     x_rtv  =  wMFISTA(y_Fista,A_Fista,At_Fista,x_rtv,lambda_k,L,mfista_iters,TV_Iters,Wx,Wy);
%     x_rtv  =  MFISTA(y_Fista,A_Fista,At_Fista,x_rtv,lambda_k,L,mfista_iters,TV_Iters);
%     
%     lambda_k = lambda_k / 2;
% %     L = L / 2;
%     
% end
% 
% 
% % x_wtv = [x_wtv(:,2:end),zeros(size(x_wtv,1),1)];
% x_rtv = [x_rtv(:,2:end),zeros(size(x_rtv,1),1)];

%% Comparing weighted to non-weighted TV
% figure;
% subplot(241);imagesc(x0);
% subplot(242);imagesc(x_tv);
% % subplot(243);imagesc(x_wtv);
% subplot(244);imagesc(x_rtv);
% subplot(246);imagesc(abs(x0-x_tv));title(num2str(psnr(x_tv,x0)));
% % subplot(247);imagesc(abs(x0-x_wtv));title(num2str(psnr(x_wtv,x0)));
% subplot(248);imagesc(abs(x0-x_rtv));title(num2str(psnr(x_rtv,x0)));
% suptitle('WMFISTA vs MFISTA');

%%  Solving using SFISTA
% x_tv_shift = [zeros(length(x_tv),1),x_tv(:,1:end-1)];
% [x_sfista, F_stack] = SFISTA( y_Fista, A_Fista, At_Fista, W, Wt, SFISTA_Params ,zeros(size(x0)), x0 );
x_sfista = x0;
% Offset of 1 pixel
% x_sfista =  [x_sfista(:,2:end), zeros(size(x_sfista,1), 1)];

%% PSNR Calculation

Err_SFISTA = psnr(x_sfista,x0)  %#ok<NOPTS>
Err_Radon  = psnr(x_radon,x0)    %#ok<NOPTS>
Err_PP     = psnr(x_cg,x0) %#ok<NOPTS>
Err_RTV    = psnr(x_radon_tv,x0)    %#ok<NOPTS>
Err_PP_TV = psnr(x_tv,x0)       %#ok<NOPTS>

%% Comaparing SFISTA with MFISTA
% figure('position',[50 50 1280 720]);
% ax(1) = subplot(231);imagesc(x0);
% title('The Phantom');
% ax(2) = subplot(232);imagesc(x_tv);
% title(['MFISTA-TV ',num2str(Err_PP_TV)]);
% ax(3) = subplot(233);imagesc(x_sfista);
% title(['SFISTA-Wavelet ',num2str(Err_SFISTA)]);
% ax(4) = subplot(234);imagesc(x_radon);
% title(['Fully Sampled Denoised FBP ',num2str(Err_Radon)]);
% ax(5) = subplot(235);imagesc(abs(x0-x_tv));
% title('Diff with MFISTA-TV');
% ax(6) = subplot(236);imagesc(abs(x0-x_sfista));
% title('Diff with SFISTA-Wavelet');
% linkaxes(ax);
% suptitle('CT PP Results');


%% Plotting

figure('Position',[90 90 1600 900]);
ax(1) = subplot(341);imagesc(x0);title('Ground Truth');

ax(2) = subplot(342);imagesc(x_radon);
title(['FBP Reconstruction ',num2str(psnr(x_radon,x0))]);
ax(3) = subplot(343);imagesc(x_cg);
title(['Pseudo-Polar CG Reocnstruction',num2str(psnr(x_cg,x0))]);
ax(4) = subplot(344);imagesc(x_tv);
title(['Pseudo-Polar with TV',num2str(psnr(x_tv,x0))]);
subplot(345);imagesc(y_radon);
title(['Equi-Angular Radon Transform ',num2str(Mproj),' Angles ',num2str(psnr(x_radon,x0))]);
ax2(1) = subplot(346);imagesc(y_pp0);title('Pseudo-Polar Sinogram');
ax2(2) = subplot(347);imagesc(ppSinogram);title('PP Sinogram from Real Projections');
ax2(3) = subplot(348);imagesc(real(invF1(D(F1(ppSinogram)))));title('Decimated PP Sinogram');
% ax2(3) = subplot(348);imagesc(abs(ppSinogram-real(invF1(A(x0)))));title('Sinograms Diff');
colorbar;

ax(5) = subplot(3,4,10);imagesc(abs(x_radon-x0));title('Diff Between Ground Truth and FBP');colorbar;
ax(6) = subplot(3,4,11);imagesc(abs(x_cg-x0));title('Diff Between Ground Truth and PP');colorbar;
ax(6) = subplot(3,4,12);imagesc(abs(x_tv-x0));title('Diff Between Ground Truth and PP');colorbar;
linkaxes(ax);linkaxes(ax2);

%% Plots for yonina - comparison
% figure;
% ax3(1) = subtightplot(131);imagesc(x0);axis off;
% title('Ground Truth');
% ax3(2) = subtightplot(132);imagesc(x_radon_tv);axis off;
% title(['FBP+TV Nyquist Sampling PSNR=',num2str(Err_RTV)]);
% ax3(3) = subtightplot(133);imagesc(x_tv);axis off;
% title(['PP sub-Nyquist PSNR=',num2str(Err_PP_TV)]);
% linkaxes(ax3);
% SaveFigure('Reconstruction',1280,440);

%% Bigger Detectors
% if strcmpi(DecMethod,'detector')
%     figure;
%     plot(real(invF1(D(F1(ppSinogram(120:390,250)),DecFactor,DecMethod))),'LineWidth',1);hold on;
%     plot(ppSinogram(120:390,250),'r','LineWidth',1);axis tight;
%     title(['Sinogram cross section at angle ',num2str(theta_pp(250))]);
%     SaveFigure('Detector',900,380);
% end

%% Sinogram Comparison
figure;
CompareImages(real(invF1(D(F1(ppSinogram)))),y_pp0,'Sinogram Comparison');
CompareImages(x_tv,x0,'PP to TV Comparison');

