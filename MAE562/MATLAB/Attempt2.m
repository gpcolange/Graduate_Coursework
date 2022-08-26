%% MAE 562 Research Paper - Proportional Navigation

clear
close all
warning off
clc

%%% Simulation Conditions

r0          = 250000                ;   % Initial Range [m]
rdot0       = 1000                  ;   % Closing Velocity [m/s]
thetadot0   = deg2rad(1)            ;   % Initial LOS angular rate [rad/s]
theta0      = deg2rad(1)            ;   % LOS angle [rad]
lambda      = 4                     ;   % Navigation Constant

time        = (0:.001:120)'                                             ;   % time [s]
IC          = [r0 rdot0 theta0 thetadot0]                                               ;   % Initial Condition Vector


options     = odeset('AbsTol',1e-8,'RelTol',1e-8,'Events', @BREAK)      ; % ODE45 solver options
[T,Z]       = ode45(@(t,z)PN(t,z,lambda) ,time,IC,options)              ;
%[T,Z]       = ode45(@(t,z)PN(t,z,lambda,r0,thetadot0,rdot0) ,time,IC,options)              ;


time1       = (0:.001:T(end))'                                          ;   
IC1         = [r0 rdot0 thetadot0]                                      ;
[T1,Z1]     = ode45(@(t,z)Free(t,z),time1,IC1,options)                  ;

% r           = Z(:,1);
% theta       = Z(:,2);
% rdot        = -sqrt(((-r0^2)*(thetadot0^2).*(Z(:,1)/r0).^(2*lambda -2)) + (rdot0)^2 + (r0^2 * thetadot0^2));
% thetadot    = thetadot0.*(Z(:,1)./r0).^(lambda -2);
 a           = 0;%lambda*sqrt((rdot).^2 + (r.^2 .* thetadot.^2)) .* thetadot;

r = Z(:,1);
rdot = Z(:,2);
theta = Z(:,3);
thetadot = Z(:,4);


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

%%% Proportional Navigation
function zdot  = PN(t,z,lambda)
z1          = z(1,1); % z1 = r
z2          = z(2,1); % z2 = rdot
z3          = z(3,1); % z3 = theta
z4          = z(4,1); % z4 = thetadot

% Equations of motion is first order form
zdot(1,1)   = z2;
zdot(2,1)   = (1-lambda)*z1*z4^2;
zdot(3,1)   = z4;
zdot(4,1)   = (lambda-2)*z2*z4/z1;
end

%%% Proportional Navigation
% function zdot  = PN(t,z,lambda,r0,thetadot0,rdot0)
% z1          = z(1,1); % z1 = r
% z2          = z(2,1); % z2 = theta
% 
% % Equations of motion is first order form
% zdot(1,1)   = sqrt((r0^2)*(thetadot0^2)*(z1/r0)^(2*lambda -2) + (rdot0)^2 + (r0^2 * thetadot0^2));
% zdot(2,1)   = thetadot0*(z1/r0)^(lambda -2);
% end


%%% Break Condition
function [value, isterminal, direction] = BREAK(t,Z)
value      = (Z(1) <= 0);
isterminal = 1;   % Stop the integration
direction  = 0;
end

%%% Plant Dynamics
function zdot  = Free(t,z)
z1          = z(1,1); % z1 = r
z2          = z(2,1); % z2 = rdot
z3          = z(3,1); % z3 = thetadot

% Equations of motion is first order form
zdot(1,1)   = z2;
zdot(2,1)   = z1*z3^3;
zdot(3,1)   = z2*z3/z1;
end

