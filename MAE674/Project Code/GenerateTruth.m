%% Propogate Quaternion from true angular rate
clear
close all
clc

% Time properties
dt          = .25;
t           = (0:dt:3600);

% True Body Angular Rate Relative to Inertial Frame, in Body Frame: Rotate about b2 axis by 90 deg/hr & b3 by 30 deg/hr 
w_B_BI      = [zeros(1,length(t)); 90*pi/(3600*180)*ones(1,length(t)); 30*pi/(3600*180)*ones(1,length(t))];

% Initialize Quaternion
q0          = (sqrt(2)/2)*[0.02,.98 , .98, 0.02]';
q0          = q0/norm(q0);
q           = zeros(4,length(t));
q(:,1)      = q0;
Mag         = zeros(1,length(t));
Mag(1)      = norm(q0);

% Numerically Intergrate Quaternion
for i = 1:length(t)-1
    q(:,i+1)    = RungeKutta4(@QuaternionKinematics,t(i),q(:,i),dt,w_B_BI(:,i));
    q(:,i+1)    = q(:,i+1)/norm(q(:,i+1));
    Mag(i+1)    = norm(q(:,i+1));
end


figure
subplot(2,1,1)
plot(t,q)
xlabel('Time [s]')
ylabel('Quaternion')
legend('q1','q2','q3','q4')
title('True Quaternion Trajectory')
grid minor


subplot(2,1,2)
plot(t,w_B_BI)
xlabel('Time [s]')
ylabel('Angular Rate [rad/s]')
legend('$\omega_1$','$\omega_2$','$\omega_3$','Interpreter','latex')
title('True $\omega_{B/I}$','Interpreter','latex')
grid minor


