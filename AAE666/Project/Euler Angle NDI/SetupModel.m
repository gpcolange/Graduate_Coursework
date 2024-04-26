clear
%% Design Contoller

% Error Dynamics: edot = A*e + B(-v), v = Ke
A               = [0 1; 0 0];
B               = [0;1];

% LQR Weighting Matrices
R               = 10;
Q               = eye(2);

% LQR Gain
K               = lqr(A,B,Q,R,[]);

use_error_loop  = 0;

%% Plant Parameters

% Moment of Inertia Matrix
Ixx             = 8800; % kg-m^2
Iyy             = 35250; % kg-m^2
Izz             = 42000; % kg-m^2
Ixz             = 2900;  % kg-m^2
I               = [Ixx 0 -Ixz;0 Iyy 0; -Ixz 0 Izz];

% Disturbance Amplitude
dist            = 0;

% Initial Conditions
psi0            = 0;
theta0          = 0;
phi0            = 0;
x0              = [phi0;theta0;psi0;0;0;0];

%% Model Parameters
Im              = I;

%% Command Generator Parameters

% Second order filter parameters
zeta            = 0.9;
wn              = pi/2;

% Create filter state space model Y/U = wn^2/(s^2 + 2*zeta*wn*s +wn^2)
A_filt          = [0 1; -wn^2 -2*zeta*wn];
B_filt          = [0;wn^2];
C_filt          = [eye(2);-wn^2 -2*zeta*wn];
D_filt          = [0;0;wn^2];

% 321 Euler Angle Commands
phi             = 30*pi/180;
psi             = 45*pi/180;
theta           = -60*pi/180;

phi_off         = 0*pi/180;
theta_off       = 0*pi/180;
psi_off         = 0*pi/180; 
