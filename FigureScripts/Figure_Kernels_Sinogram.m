%% Initializing
Initialize;

global FIGURE_PATH SaveFlag;
FIGURE_PATH = './Results/';
SaveFlag = 0;


%% Figure parameters
M = 480;
N = 480;

% Mp = 2*N+2;
Mp = N;
R = 1/pi/sqrt(2);
W = pi*M;
B = floor(Mp/8);

WindowSize = 6;
%%
t_vec     = linspace(-1/2,1/2,M);%-1/(2*size(y_rd_pad,1));
theta_vec = linspace(-1/2,1/2,N)*128;
H_full    = SinogramKernel(deg2rad(theta_vec),t_vec,B,R,W);
Window    = HammingWindow(t_vec,theta_vec/180,2*WindowSize/M);
H_window  = H_full.*Window;

H_full = H_full ./ max(H_full(:));
H_window = H_window ./ max(H_window(:));


%% The sinogram
x0   = LoadPhantom(M,'zubal');
y_rd = radon(x0,(-1/4:1/M:3/4-1/M)*180);


%%
% CompareImages(H_full,H_window);
gap = [0.0 0.008];
dr = 36;
r = N/2-dr:N/2+dr;
FigX = 1650;
FigY = FigX/6;
TitleParams = {'FontSize',34,'FontWeight','bold','FontName','Times New Roman', ...
    'Color',[1 1 0],'Units','normalized'};
label_shift = 70;
marg_h = 0;
marg_w = 0;
PosX = .025;
PosY = .94;

myfft = @(x) mat2gray(abs(fftshift(fft2(x))));

figure;
subtightplot(1,6,1,gap,marg_h,marg_w);imshow(abs(H_full(r,r)));
text(PosX,PosY,'(a)',TitleParams{:});
colormap('bone'); freezeColors;


subtightplot(1,6,3,gap,marg_h,marg_w);imshow(abs(H_window(r,r)));
text(PosX,PosY,'(c)',TitleParams{:});
colormap('bone'); freezeColors;

subtightplot(1,6,2,gap,marg_h,marg_w);imshow(myfft(H_full));
text(PosX,PosY,'(b)',TitleParams{:});
colormap('pink'); freezeColors;


subtightplot(1,6,4,gap,marg_h,marg_w);imshow(myfft(H_window));
text(PosX,PosY,'(d)',TitleParams{:});
colormap('pink'); freezeColors;

subtightplot(1,6,5,gap,marg_h,marg_w);imshow(trimcols(y_rd,N),[])
text(PosX,PosY,'(e)',TitleParams{:});
colormap('bone'); freezeColors;

subtightplot(1,6,6,gap,marg_h,marg_w);
imshow(trimcols(db(abs(fftshift(fft2(y_rd)))+1000),N),[]);colormap('pink');

colormap('pink'); freezeColors;
text(PosX,PosY,'(f)',TitleParams{:});

SaveFigure('Kernels',FigX,FigY,'notight');

%% Another figure for presentation

figure;
% subtightplot(1,6,1,gap,marg_h,marg_w);imshow(abs(H_full(r,r)));
% % text(PosX,PosY,'(a)',TitleParams{:});
% colormap('bone'); freezeColors;
% 

subtightplot(1,4,1,gap,marg_h,marg_w);imshow(abs(H_window(r,r)));
% text(PosX,PosY,'(c)',TitleParams{:});
colormap('bone'); freezeColors;

% subtightplot(1,4,1,gap,marg_h,marg_w);imshow(myfft(H_full));
% % text(PosX,PosY,'(b)',TitleParams{:});
% colormap('pink'); freezeColors;


subtightplot(1,4,2,gap,marg_h,marg_w);imshow(myfft(H_window));
% text(PosX,PosY,'(d)',TitleParams{:});
colormap('pink'); freezeColors;

subtightplot(1,4,3,gap,marg_h,marg_w);imshow(trimcols(y_rd,N),[])
% text(PosX,PosY,'(e)',TitleParams{:});
colormap('bone'); freezeColors;

subtightplot(1,4,4,gap,marg_h,marg_w);
imshow(trimcols(db(abs(fftshift(fft2(y_rd)))+1000),N),[]);colormap('pink');

colormap('pink'); freezeColors;
% text(PosX,PosY,'(f)',TitleParams{:});

SaveFigure('Kernels',FigX,FigY,'notight');