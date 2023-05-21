%% Numerical Intergrator check for van der pol oscillator

clear
close all
clc

% Model Parameters
m           = 1;
c           = 1.5;
k           = 1.2;

% time
dt          = .01;
t           = (0:dt:10);

% Initial Conditions
IC          = [1 0]';    

% ODE45 solver options
options     = odeset('AbsTol',1e-8,'RelTol',1e-8) ;

% Simulation
[T,Z]       = ode45(@(t,z) VanDerPol(t,z),t,IC,options); 

% Initialize Truth
Y           = zeros(2,length(t)) ;    
Y(:,1)      = IC;

for i = 1:length(t)-1
    Y(:,i+1)    = RungeKutta4(@VanDerPol,t(i),Y(:,i),dt);
end

figure
plot(T,Z(:,1),'-k',t,Y(1,:),'--r')
legend('ODE45','RK4')

figure
plot(T,Z(:,2),'-k',t,Y(2,:),'--r')
legend('ODE45','RK4')

%% Van Der Pol Oscillator Function
function zdot  = VanDerPol(t,z)
% Parameters
m           = 1;
c           = 1.5;
k           = 1.2;

% Define States
z1          = z(1,1); % z1 = x
z2          = z(2,1); % z2 = xdot

% Equations of motion is first order form
zdot(1,1)   = z2;
zdot(2,1)   = -2*(c/m)*(z1^2 - 1)*z2 - (k/m)*z1;
end

