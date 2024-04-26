clear
close all
clc

%% Run Model Error Loop off, perfect model
SetupModel;
simout          = sim('NDI_AttitudeControl.slx');
extract_logsout(simout)
title_str       = 'Perfect Model, No Disturbances, Controller Off';
DoPlots

%% Imperfect model with disturbance - Controller Off
psi0            = 15*pi/180;
theta0          = 20*pi/180;
phi0            = -25*pi/180;
x0              = [phi0;theta0;psi0;0;0;0];
dist            = 300;          
Im              = I + diag(.15*diag(I));
simout          = sim('NDI_AttitudeControl.slx');
extract_logsout(simout)
title_str       = 'Imperfect Model with Disturbances, Controller Off';
DoPlots

%% Imperfect model with disturbance - Controller On
use_error_loop  = 1;
simout          = sim('NDI_AttitudeControl.slx');
extract_logsout(simout)
title_str       = 'Imperfect model with Disturbances, Controller On';
DoPlots

%% Verify Regulator - Perfect model without disturbance - Controller off
% dist            = 0;    
% use_error_loop  = 0;
% Im              = I;
% phi0            = 50*pi/180;
% psi0            = 70*pi/180;
% theta0          = 30*pi/180;
% x0              = [phi0;theta0;psi0;0;0;0];
% phi             = 0;
% psi             = 0;
% theta           = 0;
% simout          = sim('NDI_AttitudeControl.slx');
% extract_logsout(simout)
% title_str       = 'Perfect model, No Disturbances, Controller Off';
% DoPlots

%% Run Model Error Loop off, perfect model with initial condition error
SetupModel;
phi_off         = 4*pi/180;
theta_off       = 8*pi/180;
psi_off         = 12*pi/180;    
x0              = [phi0 + phi_off;theta0 + theta;psi0 + psi_off;0;0;0];
simout          = sim('NDI_AttitudeControl.slx');
extract_logsout(simout)
title_str       = 'Perfect Model, with Initial Condition Error, Controller Off';
DoPlots

%% Run Model Error Loop On, perfect model with initial condition error
use_error_loop  = 1;  
simout          = sim('NDI_AttitudeControl.slx');
extract_logsout(simout)
title_str       = 'Perfect Model, with Initial Condition Error, Controller On';
DoPlots