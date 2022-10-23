%% Gabriel Colangelo MAE543 Virtual Lab1 - Liquid Level Systems
clear
close all
clc

%% Part B

% Time Vector, 10 seconds
tvec                = (0:.01:100)';
% Unit Step Input 
u                   = ones(length(tvec),1);

% Ode45 call for height
[t,height]          = ode45(@tank_control,tvec,[0],[],u,tvec,'height');

% Plot results
figure
plot(t,height)
xlabel('Time [s]')
ylabel('Height [m]')
grid minor
title('Head Height in Response to Unit Step Input')

% From Final Value Theorem for Unit Step Input f(infinity) = R [s/m^2]
R                   = height(end);

% Get Indices of Time Constant
ind_tau             = find(height >= 0.632*height(end),1);
% Time Costant
tau                 = t(ind_tau);

% From Calculated Time Constant, tau = RC
C                   = tau/R;

% Check with step function
s                   = tf('s');
sys                 = R/(R*C*s + 1);
opt                 = stepDataOptions;
opt.StepAmplitude   = 1;
y1                  = step(sys, tvec, opt);

%% Part C

% Sim time
t                   = (0:.01:2)';
% nonlinear input signal
u                   = exp(-t.^2).*cos(3.*t);

% Ode45 call for height
[Th,height]         = ode45(@tank_control,t,[0],[],u,t,'height');
% Ode45 call for flow
[Tq,flow]           = ode45(@tank_control,t,[0],[],u,t,'flow');

% H(s)/Qi(s)
H_Qi                = R/(R*C*s + 1);
% Qo(s)/Qi(s)      
Qo_Qi               = 1/(R*C*s + 1);

% lsim call for H/Qi
[yh, th]            = lsim(H_Qi,u,t);
% lsim call for Qo/Qi
[yq, tq]            = lsim(Qo_Qi,u,t);

% Plot Results
figure
plot(Th,height,'-r',th,yh,'--k')
xlabel('Time [s]')
ylabel('Head Height [m]')
grid minor
legend('ODE45','System ID')
title('$\frac{H(s)}{Q_i(s)}$ for u $= \exp(-t^2)cos(3t)$  ','interpreter','latex')

figure
plot(Tq,flow,'-r',tq,yq,'--k')
xlabel('Time [s]')
ylabel('Outflow [$\frac{m^3}{s}$]','interpreter','latex')
grid minor
legend('ODE45','System ID')
title('$\frac{Q_o(s)}{Q_i(s)}$ for u $= \exp(-t^2)cos(3t)$  ','interpreter','latex')