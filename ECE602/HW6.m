clear
close all
clc

%% Run and Plot Simulation
time        = (0:.005:1)';                              % Time [s]
IC          = [10; 10];                                 % 10 unit initial displacment at each end

% ODE45 solver options
options     = odeset('AbsTol',1e-8,'RelTol',1e-8);

% ODE45 Function call
[T, X]      = ode45(@(t,x) SuspensionSystem(t,x),time,IC,options);

% Get control input at each time
U           = zeros(size(T));
for i = 1:length(T)
    [~,U(i)]    = SuspensionSystem(T(i),X);
end

figure
subplot(3,1,1)
plot(T,X(:,1))
title('x_1 vs time')
grid minor
ylabel('x_1 displacment')
subplot(3,1,2)
plot(T,X(:,2))
title('x_2 vs time')
grid minor
ylabel('x_2 displacment')
subplot(3,1,3)
plot(T,U)
title('control input vs time')
grid minor
ylabel('u(t)')
xlabel('Time [s]')


%% State Space Function for Simple Suspension System
function [xdot, u]  = SuspensionSystem(t,x)

% State Space Matrices
A               = [-1 0; 0 -2];
B               = [1; 3];

% State Vector
x               = [x(1,1);x(2,1)]; 

% Control Law
if t <= 1
    u               = 167.83*exp(2*t - 2) - 131.47*exp(t - 1);
else 
    u = 0;
end

% State Space Equation
xdot            = A*x + B*u;
end