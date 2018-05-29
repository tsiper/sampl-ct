addpath(genpath(pwd));
close all; clear; clc;

global FIGURE_PATH SaveFlag;
FIGURE_PATH = './Results/';
SaveFlag = 0;

% Paramters
d_theta = 22.5;
d_s = 1/4;
r_points = 5;

N = 256;
Contrast = 20;

%%
theta = linspace(-45,135,2*N+2);

f = phantom(N);
p = radon(f,theta);
p_hat = real(invF1(App(f)));
p_hat = trimcols(p_hat,size(p,1));

% Resizing the Sinograms so they are square
p     = imresize(p,[2*N,2*N]);
p_hat = imresize(p_hat,[2*N,2*N]);

F = db(abs(F2(f))+Contrast);


%% Building the grids

% The polar grid
r = (-r_points:r_points)*N/(2*r_points+1);
t = -90:d_theta:90-d_theta;
[rr,tt] = meshgrid(r,t);
xx_polar = vec(rr.*cosd(tt) + N/2);
yy_polar = vec(rr.*sind(tt) + N/2);

% The psuedo polar grid
x = (-r_points+1:r_points-1)*N/(2*r_points+1);
s = -1:d_s:1;
% t_pp = [atan(s),fliplr(pi/2-atan(s))]*180/pi;
t_pp = atan(s);
% t_pp = tan(-1:d_s:1-d_s);
[rr,ss] = meshgrid(x,s);
xx_pp1 = rr;
yy_pp1 = rr.*ss;

xx_pp = vec([xx_pp1,yy_pp1])+(N)/2;
yy_pp = vec([yy_pp1,xx_pp1])+(N)/2;




%% Plotting
% Plot Options

g = 0.15;
d = 0.01;
d_text = 0.06;
PosX = .02;
PosY = .94;
% 
FigX = 1500;
FigY = 740;
% % FigX = 640;
% % FigY = 200;

scatter_opts = {'LineWidth',2,'SizeData',50,'MarkerEdgeColor',[0 0 0]};
pp_scatter = {'MarkerFaceColor',[1 1 1]};
polar_scatter = {'MarkerFaceColor',[1 1 0]};
arrow_params = {'LineWidth',8,'HeadWidth',40,'HeadStyle','vback2','HeadLength',30};
text_params  =     {'LineStyle','none',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','middle',...
    'FontWeight','bold',...
    'FontSize',36};
label_params = {'FontSize',36,'FontWeight','bold','FontName','Times New Roman', ...
    'Color',[1 1 0],'Units','normalized'};

% Plotting
h=figure('Position',[50 50 FigX FigY]);

subtightplot(2,3,[1,4],g);imshow(f,[]);
text(PosX,PosY,'(a)',label_params{:});

colormap('bone'); freezeColors;
subtightplot(2,3,2,g);imshow(p,[]);
text(PosX,PosY,'(b)',label_params{:});

colormap('bone'); freezeColors;
subtightplot(2,3,5,g);imshow(p_hat,[]);
text(PosX,PosY,'(d)',label_params{:});
colormap('bone'); freezeColors;

subtightplot(2,3,3,g);imshow(F,[]);
text(PosX,PosY,'(c)',label_params{:});
colormap('default');freezeColors;
hold on; 
scatter(xx_polar,yy_polar,scatter_opts{:},polar_scatter{:});


subtightplot(2,3,6,g);imshow(F,[]);
colormap('default');freezeColors;
text(PosX,PosY,'(e)',label_params{:});
hold on;
scatter(xx_pp,yy_pp,scatter_opts{:},pp_scatter{:});


% %% Drawing the text and arrows
% Right arrows
annotation(h,'arrow',[(2-g)/3+d (2+2*g)/3-d],[(3+g)/4 (3+g)/4],arrow_params{:});
annotation(h,'arrow',[(2-g)/3+d (2+2*g)/3-d],[(1-g)/4 (1-g)/4],arrow_params{:});

% Left Arrows
annotation(h,'arrow',[(1-2*g)/3+d (1+g)/3-d],[3/4*(1-g)/2 1/2*(1-g)/2],arrow_params{:});
annotation(h,'arrow',[(1-2*g)/3+d (1+g)/3-d],[(5+3*g)/8 (3+g)/4],arrow_params{:});

% Left Text
% Create textbox
annotation(h,'textbox',[(4+g)/6 (1-g)/4-d_text 0 0],...
    'String',{'F_1'},...
    'FontName','Mathematica5',text_params{:});

% Create textbox
annotation(h,'textbox',[(4+g)/6 (3+g)/4+d_text 0 0],...
    'String',{'F_1'},...
    'FontName','Mathematica5',text_params{:});

% Create textbox
annotation(h,'textbox',...
    [(2-g)/6 (1-g)/4 0 0],...
    'String',{'R_{pp}'},...
    text_params{:},'FontName','Mathematica5');

% Create textbox
annotation(h,'textbox',...
    [(2-g)/6 (g+3)/4 0 0],...
    'String',{'R'},...
    text_params{:},'FontName','Mathematica5');



%%% Spatial Resmpaling %%%
% Create textbox
annotation(h,'textbox',[0.5+d 0.5 0 0],...
    'String',{'Spatial','Resampling'},...
    text_params{:},'FontName','Times New Roman');
% Create arrow
annotation(h,'arrow',[(g+1)/3+d (g+1)/3+d],...
    [(1+g)/2-d (1-g)/2+d],arrow_params{:});



%%% Frequency Resampling %%%
% annotation(h,'textbox',[(5+2*g)/6+d 0.5 0 0],...
%     'String',{'Frequency','Resampling'},...
%     text_params{:},'FontName','Times New Roman');
% Create arrow
% annotation(h,'arrow',[2/3*(1+g)+d 2/3*(1+g)+d],...
%     [(1+g)/2-d (1-g)/2+d],arrow_params{:});

%% Saving the figure
SaveFigure('Transforms',FigX,FigY,'notight');
