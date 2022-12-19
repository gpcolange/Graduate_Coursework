%% Graduate Lab 4 Gabriel Colangelo
clear
close all
clc

%% Setup and Run Loop Stepper Motor Transfer Function

% Parameters
Beq         = 4e-3; 
ng          = 0.9;
nm          = 0.69;
Jeq         = 2e-3;
Jm          = 3.87e-7;
Kg          = 70;
Km          = 7.67e-3;
Kt          = 7.67e-3;
Rm          = 2.6;

% Plant transfer function theta_l/Vm = b/(a2*s^2 + a1*s + a0) 
b           = ng*Kg*nm*Kt/Rm;
a2          = Jeq;
a1          = Beq + (ng* Kg^2 * nm* Kt * Km/Rm);
a0          = 0;

% Maximum Overshoot Requirment 
Mp          = .05;

% Damping Ratio minimum
zeta_min    = sqrt(log(Mp)^2 / (log(Mp)^2 +pi^2));

% Damping ratio for minimum settling time 2% criteria (page 233)
zeta        = 0.76;

% Rise Time Requriment [s]
Tp          = .1;

% Natural frequency [rad/s]
wn          = pi/(Tp*sqrt(1 - zeta^2));

% Closed form solution for PD control gains 
Kp          = (a2/b)*wn^2;              % [V/rad]
Kd          = (2*zeta*wn*a2 - a1)/b;    % [V/rad/s]

time        = 0:.001:5;  % Sim Time Vector
output      = sim('StepperMotor.slx',time);

overshoot   = output.LoadPositionStep(end) + Mp;

figure
plot(output.tout,output.LoadPositionStep,output.tout,output.StepCommand)
xline(1+Tp,'--k',{'Peak','Time','Requirement'});
yline(overshoot,'--r',{'Overshoot','Requirement'});
legend('Response','Command','Location','southeast')
ylabel('Load Position [rad]')
grid minor
xlabel('Time [s]')
titlestr    = ['Step Response using Closed Loop Poles: $\zeta$ = ', num2str(zeta), ', $K_p$ = ', num2str(Kp),...
               ' [V/rad], $K_d$ = ',num2str(Kd),' [V/rad-s]'];
title(titlestr,'Interpreter','latex')

%% Tune Controller
Kp          = Kp - 8.5;

time        = 0:.001:5;  % Sim Time Vector
output      = sim('StepperMotor.slx',time);

figure
plot(output.tout,output.LoadPositionStep1,output.tout,output.StepCommand1)
xline(1+Tp,'--k',{'Peak','Time','Requirement'});
yline(overshoot,'--r',{'Overshoot','Requirement'});
legend('Response','Command','Location','southeast')
ylabel('Load Position [rad]')
grid minor
xlabel('Time [s]')
titlestr    = ['Step Response with Tuned Control Gains: $K_p$ = ', num2str(Kp),...
               ' [V/rad], $K_d$ = ',num2str(Kd),' [V/rad-s]'];
title(titlestr,'Interpreter','latex')


% Verify Tuned Controller Meets Specs
s           = tf('s');
G_CL        = b*(Kp + Kd*s)/(a2*s^2 + (a1 + b*Kd)*s + b*Kp);

ind_max     = find(output.LoadPositionStep == max(output.LoadPositionStep),1);
Mp_calc     = (output.LoadPositionStep(ind_max) - output.LoadPositionStep(end))/output.LoadPositionStep(end) *100
tp_calc     = output.tout(ind_max) - 1

figure
plot(output.tout,output.LoadPositionSine,output.tout,output.SineCommand)
legend('Response','Command','Location','southeast')
ylabel('Load Position [rad]')
grid minor
xlabel('Time [s]')
title('Sine Wave Input')

figure
plot(output.tout,output.LoadPositionSquare,output.tout,output.SquareCommand)
legend('Response','Command')
ylabel('Load Position [rad]')
grid minor
xlabel('Time [s]')
title('Square Wave Input')

figure
nichols(G_CL)
ngrid
