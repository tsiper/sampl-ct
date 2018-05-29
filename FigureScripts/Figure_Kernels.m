%% Initializing
Initialize;
set(0,'DefaultTextInterpreter','none');

%% Figure parameters
M = 512;
N = 512;
WideFlag = 1;

if WideFlag
    sub_a = 1; sub_b = 4;
    FigX = 1280;
    FigY = 320;
    FontSize = 48;
    FileName = 'Kernels_Wide';
else
    sub_a = 2; sub_b = 2;
    FigX = 640;
    FigY = 640;
    FontSize = 34;
    FileName = 'Kernels';
end


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


%%
% CompareImages(H_full,H_window);
gap = [0.01 0.01];
dr = 36;
r = N/2-dr:N/2+dr;

% TitleParams = {'FontSize',34,'FontWeight','bold','FontName','Times New Roman', ...
%     'Color',[1 1 0],'Units','normalized'};
TitleParams = {'FontSize',FontSize,'FontWeight','bold', ...
    'FontName','Times New Roman','Color',[1 1 0],'Units','normalized'};
label_shift = 70;
marg_h = 0;
marg_w = 0;
PosX = .025;
PosY = .9;



myfft = @(x) mat2gray(abs(fftshift(fft2(x))));

figure('Position',[50 50 FigX FigY]);
subtightplot(sub_a,sub_b,1,gap,marg_h,marg_w);imagesc(abs(H_full(r,r))); axis off;
text(PosX,PosY,'(a)',TitleParams{:});
colormap bone; freezeColors;

subtightplot(sub_a,sub_b,3,gap,marg_h,marg_w);imagesc(abs(H_window(r,r)));axis off;
text(PosX,PosY,'(c)',TitleParams{:});
colormap bone; freezeColors;

subtightplot(sub_a,sub_b,2,gap,marg_h,marg_w);imagesc(myfft(H_full));axis off;
text(PosX,PosY,'(b)',TitleParams{:});
colormap pink; freezeColors;

subtightplot(sub_a,sub_b,4,gap,marg_h,marg_w);imagesc(myfft(H_window));axis off;
text(PosX,PosY,'(d)',TitleParams{:});
colormap pink; freezeColors;


%% Saving
SaveFigure(FileName,FigX,FigY);