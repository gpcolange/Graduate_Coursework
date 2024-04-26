clear
close all
clc

%% Run Model Error Loop off, perfect model
% Setup and Run Sim
SetupModel;
simout          = sim('NDI_model.slx');
extract_logsout(simout)
title_str       = 'Perfect Model, No Disturbances, Controller Off';
DoPlots

%% Add sinusoidal disturbance torque, change IC
% Initial Conditions
psi0            = 15*pi/180;
theta0          = 20*pi/180;
phi0            = -25*pi/180;
q0              = Euler321toQuat(phi0, theta0, psi0);
x0              = [q0;0;0;0];
dist            = 500;          

simout          = sim('NDI_model.slx');
extract_logsout(simout)
title_str       = 'Perfect Model with Disturbances, Controller Off';
DoPlots

%% Imperfect model, no disturbance
% Initial Conditions
dist            = 0;          

% Create errors in inertia modeling
Im              = I + 300*randn(3,3);
simout          = sim('NDI_model.slx');
extract_logsout(simout)
title_str       = 'Imperfect model, No Disturbances, Controller Off';
DoPlots

%% Imperfect model with disturbance - Controller On
dist            = 250;    
use_error_loop  = 1;
simout          = sim('NDI_model.slx');
extract_logsout(simout)
title_str       = 'Imperfect model with Disturbances, Controller On';
DoPlots

%% Verify Regulator - Perfect model without disturbance - Controller off
dist            = 0;    
use_error_loop  = 0;
Im              = I;
phi0            = 50*pi/180;
psi0            = 70*pi/180;
theta0          = 30*pi/180;
q0              = Euler321toQuat(phi0, theta0, psi0);
x0              = [q0;0;0;0];
phi             = 0;
psi             = 0;
theta           = 0;
simout          = sim('NDI_model.slx');
extract_logsout(simout)
title_str       = 'Perfect model, No Disturbances, Controller Off';
DoPlots