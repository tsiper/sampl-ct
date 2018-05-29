%% Initializing
clear;
close all;
clc;

global FIGURE_PATH SaveFlag;
FIGURE_PATH = './Results/';
SaveFlag = 0;

addpath './Utilities';

%% Creating the grids
% the cartesian grid
x = -4:4;
% Cartesian plot
[xx,yy]=meshgrid(x,x);

% The polar grid
d_theta = 15;
r = -4:4;
t = -90:d_theta:90-d_theta;
[rr,tt] = meshgrid(r,t);
xx_polar = rr.*cosd(tt);
yy_polar = rr.*sind(tt);

% The psuedo polar grid
d_s = 1/4;
s = -1:d_s:1;
% t_pp = [atan(s),fliplr(pi/2-atan(s))]*180/pi;
t_pp = atan(s);
% t_pp = tan(-1:d_s:1-d_s);
[rr,ss] = meshgrid(x,s);
xx_pp = rr;
yy_pp = rr.*ss;

xx_pp1 = [xx_pp,yy_pp];
yy_pp1 = [yy_pp,xx_pp];



%% Plot Paramters
% the arrows
xO = 0.3;  
yO = 0.15;

LineWidth = 2;
FigX = 1340;
FigY = 420;
TitleFontSize = 24;
gap = 0.02;
AxisLineWidth = 2;

%% PLOTTING

h = figure('Color',[1 1 1],'Position',[50 50 FigX FigY]);

subtightplot(1,3,1,gap);
scatter(xx(:),yy(:),'LineWidth',LineWidth,'MarkerEdgeColor',[0.8 0 0]);
axis([-5 5 -5 5])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
set(ax, 'Xtick', [], 'Ytick', [], 'box', 'off','LineWidth',AxisLineWidth)
patch(...
    [5-xO -yO; 5-xO +yO; 5.0 0.0], ...
    [yO 5-xO; -yO 5-xO; 0 5], 'k', 'clipping', 'off')
title('Cartesian Grid','FontSize',TitleFontSize);


subtightplot(1,3,3,gap);
scatter(xx_polar(:),yy_polar(:),'LineWidth',LineWidth,'MarkerEdgeColor',[0 0 0.8]);
axis([-5 5 -5 5])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
set(ax, 'Xtick', [], 'Ytick', [], 'box', 'off','LineWidth',AxisLineWidth);
patch(...
    [5-xO -yO; 5-xO +yO; 5.0 0.0], ...
    [yO 5-xO; -yO 5-xO; 0 5], 'k', 'clipping', 'off')
title('Polar Grid','FontSize',TitleFontSize);


subtightplot(1,3,2,gap);
scatter(xx_pp1(:),yy_pp1(:),'LineWidth',LineWidth,'MarkerEdgeColor',[0 0.4 0]);
axis([-5 5 -5 5])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
set(ax, 'Xtick', [], 'Ytick', [], 'box', 'off','LineWidth',AxisLineWidth);
patch(...
    [5-xO -yO; 5-xO +yO; 5.0 0.0], ...
    [yO 5-xO; -yO 5-xO; 0 5], 'k', 'clipping', 'off')
title('Pseudo-Polar Grid','FontSize',TitleFontSize);

%% Darwing the artwork
% annotation(h,'arrow',[0.625 0.725],...
%     [0.25 0.25],'LineWidth',6);


%% Saving the figure!
SaveFigure('Grids',FigX,FigY);