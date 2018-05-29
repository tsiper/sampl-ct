%% Initializing
clear; close all; clc;
global FIGURE_PATH SaveFlag;

% Making sure we use regular font
set(0,'defaulttextinterpreter','none');

%% Paramters

% The place we save the figures to
FIGURE_PATH = './';

% Should we save the figure? Set to 1 if yes
SaveFlag = 1;

% Uncheck this to save the figures
FigureName = 'MyPhantom';

%% Generating the pictures
x = cell(1,6);
x{1} = phantom(256);
sig = 0.02;

for i=2:6
    x{i} = x{1}+i*sig*randn(size(x{1}));
end
    


%% Plotting

Width = 600; Height = 900;
figure('Position',[50 50 Width Height]);
g = 0.01; h = 0; w = 0;

ax(1) = subtightplot(3,2,1,g,h,w);imshowzoom(x{1},'a');
ax(2) = subtightplot(3,2,2,g,h,w);imshowzoom(x{2},'b');
ax(3) = subtightplot(3,2,3,g,h,w);imshowzoom(x{3},'c');
ax(4) = subtightplot(3,2,4,g,h,w);imshowzoom(x{4},'d');
ax(5) = subtightplot(3,2,5,g,h,w);imshowzoom(x{5},'e');
ax(6) = subtightplot(3,2,6,g,h,w);imshowzoom(x{6},'f');
colormap('bone');
linkaxes(ax);

SaveFigure(FigureName);

