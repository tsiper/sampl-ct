%% Initializies for a clean run

% Going to the root directory
cd(fileparts(which(mfilename)));

% Adding the entire path
addpath(genpath(pwd));

% Cleaning stuff up
close all;
clear;

%% Some default styling for publishing
set(0,'defaulttextinterpreter','tex')
set(groot, 'defaultAxesTickLabelInterpreter','tex'); 
set(groot, 'defaultLegendInterpreter','tex');

%% Setting up the globals I like to use
global SaveFlag FIGURE_PATH PlotFlag DebugFlag; 
FIGURE_PATH = './Results/'; 
SaveFlag = 0;
PlotFlag = 0;
DebugFlag = 0;
