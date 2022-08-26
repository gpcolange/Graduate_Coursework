%% MAE 562 Research Paper - Proportional Navigation

clear
close all
clc

%%% Simulation Conditions

r0          = 150000                ;   % Initial Range [m]
rdot0       = 1500                  ;   % Closing Velocity [m/s]
thetadot0   = deg2rad(1)            ;   % Initial LOS angular rate [rad/s]
theta0      = deg2rad(10)           ;   % LOS angle [rad]
lambda      = 4                     ;   % Navigation Constant

time        = (0:.001:120)'                                             ;   % time [s]
IC          = [r0 rdot0 theta0 thetadot0]                               ;   % Initial Condition Vector


options     = odeset('AbsTol',1e-8,'RelTol',1e-8,'Events', @BREAK)      ; % ODE45 solver options
[T,Z]       = ode45(@(t,z)PN(t,z,lambda) ,time,IC,options)              ;

options1    = odeset('AbsTol',1e-8,'RelTol',1e-8)                       ;
time1       = (0:.001:T(end))'                                          ;   
IC1         = [r0 rdot0 thetadot0]                                      ;
[T1,Z1]     = ode45(@(t,z)Free(t,z),time1,IC1,options1)                 ;

r           = Z(:,1);
rdot        = Z(:,2);
theta       = Z(:,3);
thetadot    = Z(:,4);

a           = max(abs(lambda*sqrt((r.*thetadot.^2).^2 + (rdot.*thetadot).^2)));
avec        = lambda*sqrt((r.*thetadot.^2).^2 + (rdot.*thetadot).^2);
v           = sqrt(rdot.^2 + r.^2.*thetadot.^2);

figure
ax1 = subplot(4,1,1);
plot(T,r)
xlabel('Time [s]')
ylabel('Distance to Target [m]')
grid minor
title({['With Proportional Navigation'],['Max Acceleration: ',num2str(max(abs(a))),' [m/s^2]']})

ax2 = subplot(4,1,2);
plot(T,rdot)
xlabel('Time [s]')
ylabel('Closing velocity [m/s]')
grid minor

ax3 = subplot(4,1,3);
plot(T,rad2deg(thetadot))
xlabel('Time [s]')
ylabel('Line of Sight Rate [deg/s]')
grid minor

ax4 = subplot(4,1,4);
plot(T,rad2deg(theta))
xlabel('Time [s]')
ylabel('Line of Sight [deg]')
grid minor
linkaxes([ax1 ax2 ax3 ax4],'x');

figure
h1 = subplot(2,1,1);
plot(T,avec)
xlabel('Time [s]')
ylabel('Command Acceleration [m/s^2]')
grid minor
title('Ideal Proportional Navigation Acceleration and Velocity')

h2 = subplot(2,1,2);
plot(T,v)
xlabel('Time [s]')
ylabel('Velocity Magnitude [m/s]')
ylim([max(v)/2, 2*max(v)])
grid minor
linkaxes([h1 h2],'x');

figure
ax1 = subplot(3,1,1);
plot(T1,Z1(:,1))
xlabel('Time [s]')
ylabel('Distance to Target [m]')
grid minor
title('No Guidance')

ax2 = subplot(3,1,2);
plot(T1,Z1(:,2))
xlabel('Time [s]')
ylabel('Closing velocity [m/s]')

grid minor
ax3 = subplot(3,1,3);
plot(T1,rad2deg(Z1(:,3)))
xlabel('Time [s]')
ylabel('Line of Sight Rate [deg/s]')
grid minor
linkaxes([ax1 ax2 ax3],'x');

%% Generate Table Values
for lambda = 3:7
    options     = odeset('AbsTol',1e-8,'RelTol',1e-8,'Events', @BREAK)      ; % ODE45 solver options
    [T,X]       = ode45(@(t,z)PN(t,z,lambda) ,time,IC,options)              ;
    
    r           = Z(:,1);
    rdot        = Z(:,2);
    theta       = Z(:,3);
    thetadot    = Z(:,4);
    a           = max(abs(lambda*sqrt((r.*thetadot.^2).^2 + (rdot.*thetadot).^2)));
    
    fprintf('PN Constant = %i \n',lambda);
    fprintf('Intercept Time: %.2f [s] \n',T(end));
    fprintf('Max acceleration = %.2f [m/s^2] \n',a);
    fprintf('\n');
end


