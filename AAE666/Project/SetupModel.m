clear
close all
clc

%% Design Contoller

% Error Dynamics: edot = A*e + B(-v), v = Ke
A               = [0 1; 0 0];
B               = [0;1];

% LQR Weighting Matrices
R               = 10;
Q               = .1*eye(2);

% LQR Gain
K               = lqr(A,B,Q,R,[]);

use_error_loop  = 0;

%% Plant Parameters

% Moment of Inertia Matrix
Ixx             = 12900; % kg-m^2
Iyy             = 75500; % kg-m^2
Izz             = 85500; % kg-m^2
Ixz             = 1300;  % kg-m^2
I               = [Ixx 0 -Ixz;0 Iyy 0; -Ixz 0 Izz];

% Disturbance Amplitude
dist            = 0;

% Initial Conditions
x0              = [0,0,0,1,0,0,0]';

%% Model Parameters
Im              = I;

%% Command Generator Parameters

% Second order filter parameters
zeta            = 0.9;
wn              = pi/2;

% Create filter state space model Y/U = wn^2/(s^2 + 2*zeta*wn*s +wn^2)
A               = [0 1; -wn^2 -2*zeta*wn];
B               = [0;wn^2];
C               = [eye(2);-wn^2 -2*zeta*wn];
D               = [0;0;wn^2];

% 321 Euler Angle Commands
phi             = 30*pi/180;
psi             = 45*pi/180;
theta           = -60*pi/180;
