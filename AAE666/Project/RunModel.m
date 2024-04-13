clear
close all
clc

%% Run Model Error Loop off, perfect model
% Setup and Run Sim
SetupModel;
simout = sim('NDI_model.slx');

% Extract data and plot
extract_logsout(simout)
title_str1 = 'Quaternions: Perfect Model, No Disturbances, Controller Off';
title_str2 = 'Euler Angles: Perfect Model, No Disturbances, Controller Off';
DoPlots

