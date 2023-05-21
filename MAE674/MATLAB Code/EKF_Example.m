%% EKF Example for Van der Pol Oscillator

clear
close all
clc

%% Generate Measurements
% time
dt          = .01;
t           = (0:dt:10);

% Initial Conditions
X           = zeros(2,length(t)) ;    
X(:,1)      = [1 0]';

for i = 1:length(t)-1
    X(:,i+1)    = RungeKutta4(@VanDerPol,t(i),X(:,i),dt);
end

% Measurement Error standard deviation
std_m       = .01;

% Position Measurements - Y = [1 0]x
Ym          = X(1,:) + std_m*randn(1,length(X));

%% Initialize 
% Initialize Error Covariance
P0          = 1000*eye(2);
Pcov        = zeros(2,length(t));
Pcov(:,1)   = diag(P0);
P           = P0;

% Initialize Estimates
xe          = zeros(2,length(t));
xe(:,1)     = [1 0]';

m           = 1;
c           = 1.5;
k           = 1.2;

for i = 1:length(t)-1 
    %% Gain
    H               = [1 0];
    h               = xe(1,i);
    % Measurment Error Spectral Density
    R               = std_m^2*eye(size(Ym(i)*Ym(i)'));
    % Kalman Gain
    K               = P*H'*inv(H*P*H' + R);

    %% Update
    xe(:,i)         = xe(:,i) + K*(Ym(i) - h);
    P               = (eye(size(P)) - K*H)*P;

    %% Propogate
    xe(:,i+1)       = RungeKutta4(@VanDerPolError,t(i),xe(:,i),dt);
    % Linearized Dynamics 
    F               = [0, 1; -4*(c/m)*xe(1,i)*xe(2,i) - (k/m), -2*(c/m)*(xe(1,i)^2 - 1)];
    G               = [0 1]';
    % Process Noise Spectral Density
    Q               = 0.2*[0 0; 0 1]; 
    % Continous to Discrete
    Phi             = expm(F*dt);

    % Use discrete time solution for error covariance
    P               = Phi*P*Phi' + (Q*dt);
    Pcov(:,i+1)     = diag(P);
end

figure
plot(t,Ym,t,X(1,:))
xlabel('Time [s]')
ylabel('Position [m]')
grid minor
legend('Measurement','Truth')

figure
plot(t,X(1,:),t,xe(1,:))
xlabel('Time [s]')
ylabel('Position [m]')
grid minor
legend('Truth','Estimate')

figure
plot(t,X(2,:),t,xe(2,:))
xlabel('Time [s]')
ylabel('Velocity [m/s]')
grid minor
legend('Truth','Estimate')

figure
plot(t,X(1,:) - xe(1,:),t,3*sqrt(Pcov(1,:)),'--r',t,-3*sqrt(Pcov(1,:)),'--r')
xlabel('Time [s]')
ylim([-.1 .1])
ylabel('Position Errors [m]')
grid minor

figure
plot(t,X(2,:) - xe(2,:),t,3*sqrt(Pcov(2,:)),'--r',t,-3*sqrt(Pcov(2,:)),'--r')
xlabel('Time [s]')
ylim([-.5 .5])
ylabel('Velocity Errors [m]')
grid minor

%% Van Der Pol Oscillator Function
function zdot  = VanDerPol(t,z)
% Parameters
m           = 1;
c           = 1;
k           = 1;

% Define States
z1          = z(1,1); % z1 = x
z2          = z(2,1); % z2 = xdot

% Equations of motion is first order form
zdot(1,1)   = z2;
zdot(2,1)   = -2*(c/m)*(z1^2 - 1)*z2 - (k/m)*z1;
end

function zdot  = VanDerPolError(t,z)
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