%% This Problem Performs Poorly For EKF
clear
close all
clc

%% Generate Measurements
dt      = 1;
tm      = 0:dt:60;
t       = 0:.01:60;
a       = 5e-5;     % Air density constant 

% Initial Conditions
X       = zeros(3,length(tm)) ;   
X(:,1)  = [3e5, 2e4, 1e-3]';

ytilde  = zeros(1,length(tm));

% Horizontal Distance
M       = 1e5;

% Vertical Distance
Z       = M;

% Variance
r       = 1e4;

for i = 1:length(tm)-1
    X(:,i+1)    = RungeKutta4(@FallingBody,tm(i),X(:,i),dt);  
end

ytilde  = sqrt(M^2 + (X(1,:) - Z).^2) + sqrt(r)*randn(size(ytilde));

xe      = zeros(size(X));
xe(:,1) = [3e5, 2e4, 3e-5]';

P0      = diag([1e6, 4e6, 1e-4]);

% Initiliaze State Error Covariance
P           = P0;
Pcov        = zeros(3,length(t));
Pcov(:,1)   = diag(P);

j           = 0;

%% Extended Kalman Filter
for i = 1:length(t)-1 
   
    %%% Update
    ind             = find(t(i) == tm,1);
    if ~isempty(ind)
        j           = j + 1;

        %%% Gain

        % H = dh/dx
        H               = [((xe(1,i) - Z)/sqrt(M^2 + (xe(1,i) - Z)^2)), 0, 0];
        h               = sqrt(M^2 + (xe(1,i) - Z)^2);
    
        % Measurment Error Spectral Density
        R               = r;
        % Kalman Gain
        K               = P*H'*inv(H*P*H' + R);
    
        xe(:,i)         = xe(:,i) + K*(ytilde(j) - h);
        P               = (eye(size(P)) - K*H)*P;
    end

    %%% Propogate
    xe(:,i+1)       = RungeKutta4(@FallingBody,t(i),xe(:,i),.01);

    % Linearized Dynamics 
    F               = exp(-a*xe(1,i))*[0,-exp(a*xe(1,i)), 0;...
                      a*xe(2,i)^2*xe(3,i), -2*xe(2,i)*xe(3,i), -xe(2,i)^2;...
                      0, 0, 0];

    % Process Noise Covariance
    Qk              = 0;
    
    % Continous to Discrete
    Phi             = expm(F*dt);

    % Use discrete time solution for error covariance
    P               = Phi*P*Phi' + Qk;
    Pcov(:,i+1)     = diag(P);
end

figure
subplot(3,1,1)
plot(tm,X(1,:),t,xe(1,:))
xlabel('Time [s]')
ylabel('Altitude [m]')
grid minor
legend('Truth','Estimate')
title('EKF Altitude Estimate')
subplot(3,1,2)
plot(tm,X(2,:),t,xe(2,:))
xlabel('Time [s]')
ylabel('Velocity [m/s]')
grid minor
title('EKF Velocity Estimate')
subplot(3,1,3)
plot(tm,X(3,:),t,xe(3,:))
xlabel('Time [s]')
ylabel('Ballisitic Coefficient')
grid minor
title('EKF Ballisitic Coefficient Estimate')

function xdot = FallingBody(t,x)
    a           = 5e-5;     % Air density constant 
    
    x1          = x(1,1);   % Altitude State
    x2          = x(2,1);   % Downward Velocity State
    x3          = x(3,1);   % Ballistic Coefficient
    
    xdot(1,1)   = -x2;
    xdot(2,1)   = -exp(-a*x1)*x2^2*x3;
    xdot(3,1)   = 0;
end