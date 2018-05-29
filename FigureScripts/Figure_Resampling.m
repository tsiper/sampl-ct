%% Resampling Radon to PP using our specific sub-space Kernel

%% Initializing
Initialize;

SaveFlag = 0;
PlotFlag = 1;
DebugFlag = 1;

%% PARAMETERS
N  = 256;                  % The image side resolution

% Chossing a phantom
PhantomTypes = {'brain','thorax','shepp','zubal'};
PhantomType  = PhantomTypes{4};
ScanType     = 'MATLAB';

WindowSize =  6;       % The radius of the window
Sigma      = 0.015;     % Noise variance
Mp         = 2*N + 2;  % Number of sensors
cg_iters   = 3;        % Maximum iterations for conjugate gradient alg.
DecFactor  = 1;        % Decimation factor of the measurements

% The default parameters for the kernel
B         = floor(Mp/10);
% R         = 1/pi/2/sqrt(2);
R         = 1/pi/2/sqrt(2);
W         = (pi)*(2*N+1);

%% Generating the scanned data
% Loading the Phantom
x0 = LoadPhantom(N,PhantomType);

% Number of sensors
N_rd     = 2*ceil(sqrt(2)*(N-floor((N-1)/2)-1))+3;

% Designing the scanning angles
theta      = (-1/4:1/Mp:3/4-1/Mp)*180;
theta_pad  = (-1/2:1/Mp:1-1/Mp)*180;
% theta_pad = theta;

% The desired PP scanning angles
l = -N/2:N/2;
theta_pp = [atan(2*l/N),fliplr(pi/2-atan(2*l/N))]*180/pi;

% Scanning the object using regular Radon transform
if strcmpi(ScanType,'AIR')
    A_rd     = paralleltomo(N,theta,N_rd);
    A_rd_pad = paralleltomo(N,theta_pad,N_rd);
    y_rd     = reshape(A_rd*x0(:),[N_rd,length(theta)]);
    y_rd_pad     = reshape(A_rd*x0(:),[N_rd,length(theta)]);
elseif strcmpi(ScanType,'MATLAB')
    y_rd     = radon(x0,theta);
    y_rd_pad = radon(x0,theta_pad);
    y_rd_pad = padcols(y_rd_pad,2*N+1);
else
    error('Please choose ''ScanType'' as either ''AIR'' or ''MATLAB''');
end
    
y_rd_pad   = y_rd_pad+Sigma*range(y_rd_pad(:))*randn(size(y_rd_pad));
y_rd_pad(y_rd_pad<0) = 0;
y_rd_noise = y_rd +Sigma*(max(y_rd(:))-min(y_rd(:)))*randn(size(y_rd));
y_rd_noise(y_rd_noise<0) = 0;
y_pp0 = real(invF1(App(x0)));

%% Decimating y_rd_pad
y_rd_pad  = y_rd_pad(:,1:DecFactor:end);
theta_pad = theta_pad(1:DecFactor:end);

%% PP Interpolation
[ y_pp_filt ]   = InterpolateWindow( y_rd_pad, N, theta_pad , WindowSize, B,R,W);
[ y_pp_spline ] = InterpolateSinogram( y_rd_pad, N, theta_pad, 'spline' );
[ y_pp_linear ] = InterpolateSinogram( y_rd_pad, N, theta_pad,'linear' );


% Moving to the positive side
y_pp_filt(y_pp_filt<0) = 0;
y_pp_spline(y_pp_spline<0) = 0;
y_pp_linear(y_pp_linear<0) = 0;

% Normalizing
y_pp_filt   = y_pp_filt   * fminsearch(@(Amp) -psnr(Amp*(y_pp_filt),y_pp0),1);
y_pp_spline = y_pp_spline * fminsearch(@(A) -psnr(A*(y_pp_spline),y_pp0),1);
y_pp_linear = y_pp_linear * fminsearch(@(A) -psnr(A*(y_pp_linear),y_pp0),1);

%% Reconstructions
x_pp     = Proj_C(conjgrad(@(x) App_T(M(App(x))), App_T(M(F1(y_pp_filt))) , zeros(N) ,cg_iters ));
x_spline = Proj_C(conjgrad(@(x) App_T(M(App(x))), App_T(M(F1(y_pp_spline))) , zeros(N) ,cg_iters ));
x_linear = Proj_C(conjgrad(@(x) App_T(M(App(x))), App_T(M(F1(y_pp_linear))) , zeros(N) ,cg_iters ));

if strcmpi(ScanType,'MATLAB')
    x_fbp    = iradon(y_rd_noise,theta);
    x_fbp    = x_fbp(2:end-1,2:end-1);
else
    x_fbp    = Proj_C(vec2im(fbp(A_rd,y_rd_noise,theta)));
end

x_pp(x_pp<0) = 0;
x_spline(x_spline<0) = 0;
x_linear(x_spline<0) = 0;

x_pp = x_pp * fminsearch(@(A) -psnr(A*x_pp,x0),1);
x_spline = x_spline * fminsearch(@(A) -psnr(A*x_spline,x0),1);
x_linear = x_linear * fminsearch(@(A) -psnr(A*x_linear,x0),1);

%% Calculations
yPSNR_spline = psnr(y_pp_spline,y_pp0);
yPSNR_linear = psnr(y_pp_linear,y_pp0);
yPSNR_filt   = psnr(y_pp_filt,y_pp0);
yPSNR_fbp    = psnr(y_rd_noise,y_rd);

xPSNR_spline = psnr(x_spline,x0);
xPSNR_linear = psnr(x_linear,x0);
xPSNR_filt   = psnr(x_pp,x0);
xPSNR_fbp    = psnr2(x_fbp,x0);


%% Plotting
CompareImages(x_spline,x0, 'Spline Recon.');
CompareImages(x_fbp,x0,    'FBP Recon.');
CompareImages(x_pp,x0,   'Filt Recon.');
CompareImages(x_linear,x0,'Linear Recon.');
%% Sinogram Plotting
CompareImages(y_pp_spline,y_pp0,  'Spline Sinogram');
CompareImages(y_pp_filt,y_pp0,    'Filt. Sinogram');

y_pp0_trim       = trimcols(y_pp0,N_rd);
y_pp_filt_trim   = trimcols(y_pp_filt,N_rd);
y_pp_spline_trim = trimcols(y_pp_spline,N_rd);

%% Figure 1 - Showing the polar sinograms
gap = [.02 .01];
fft_contrast = 3000;
marg_h = 0;
marg_w = 0;
FigX = 1500;
FigY = 500;
TitleParams = {'FontSize',60,'FontWeight','bold','FontName','Times New Roman', ...
    'Color',[1 1 0],'Units','normalized'};
PosX = .025;
PosY = .92;

ColormapSino  = 'bone'; ColormapRecon = 'bone'; ColormapFreq  = 'pink';

myimshow = @(x) imagesc(mat2gray(x));
myfft    = @(x) mat2gray(db(abs(fftshift(fft2(x)))+fft_contrast));

figure('Name','Polar','Position',[50 50 FigX FigY],'units','normalized');

subtightplot(1,3,1,gap,marg_h,marg_w);  myimshow(x0);
axis off square; text(PosX,PosY,'(a)',TitleParams{:});
colormap(ColormapSino); freezeColors;
% subtightplot(2,3,4,gap,marg_h,marg_w); myimshow(y_rd_noise);
% text(PosX,PosY,'(d)',TitleParams{:}); axis off;
% colormap(ColormapSino); freezeColors;
subtightplot(1,3,2,gap,marg_h,marg_w); myimshow(y_rd);
text(PosX,PosY,'(b)',TitleParams{:}); axis square off;
colormap(ColormapRecon); freezeColors;
% subtightplot(2,3,5,gap,marg_h,marg_w); myimshow(myfft(y_rd_noise));
% text(PosX,PosY,'(e)',TitleParams{:}); axis off;
% colormap(ColormapFreq); freezeColors;
subtightplot(1,3,3,gap,marg_h,marg_w);myimshow(myfft(y_rd));
text(PosX,PosY,'(c)',TitleParams{:}); axis square off;
colormap(ColormapFreq); freezeColors; 
% subtightplot(2,3,6,gap,marg_h,marg_w); myimshow(x_fbp);
% text(PosX,PosY,'(f)',TitleParams{:}); axis off;
% colormap(ColormapRecon); freezeColors;

SaveFigure('Polar',FigX,FigY,'notight');

%% Figure 2 - Showing the PP results
figure('Position',[50 50 FigX 3/2*FigY]);

colormap('bone');

subtightplot(231,gap); myimshow(y_pp0_trim);        axis off;
text(PosX,PosY,'(a)',TitleParams{:});

subtightplot(232,gap); myimshow(y_pp_spline_trim);  axis off;
text(PosX,PosY,'(b)',TitleParams{:});

subtightplot(233,gap); myimshow(y_pp_filt_trim);    axis off;
text(PosX,PosY,'(c)',TitleParams{:});

% subtightplot(234,gap); myimshow(x0);                axis off;
% text(PosX,PosY,'(d)',TitleParams{:});
% 
subtightplot(235,gap); imagesc(abs(y_pp0_trim-y_pp_spline_trim),[0 10]);          axis off;
text(PosX,PosY,'(e)',TitleParams{:});
colorbar;

subtightplot(236,gap); imagesc(abs(y_pp0_trim-y_pp_filt_trim),[0 10]);          axis off;
text(PosX,PosY,'(e)',TitleParams{:});
colorbar;

% SaveFigure('PseudoPolar',FigX,FigY,'notight');


