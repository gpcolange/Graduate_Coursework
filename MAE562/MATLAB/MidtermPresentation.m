%% Midterm Presentation Script
clear
close all 
clc

%% Reachable Set

%%% Target Trajectory

v           = -2058;       % [m/s] Target Velocity - Mach 6
h0          = 10000;       % [m] Target Initial Height
x0          = 100000;      % [m] Target horizontal distance
theta       = 2;           % [deg] Flight path angle

t           = (0:.001:133)'; % time [s]


y           = v*sind(theta)*t + h0;
x           = v*cosd(theta)*t + x0;

%%% Neccessary Interceptor Velocity

intercept_height    = 500; % [m]
ind                 = find(y < intercept_height,1);
            
theta_I             = atan2d(y(ind),x(ind));
vI                  = x(ind)/(t(ind)*cosd(theta_I));

yI                  = vI*sind(theta_I)*t;
xI                  = vI*cosd(theta_I)*t;

plot(t,sqrt((x - xI).^2 + (y - yI).^2))
xlabel('Time [s]')
ylabel('Distance [m]')
title('Difference in Position')
grid minor

figure
plot3(t,xI,yI,t,x,y)
xlabel('Time [s]')
ylabel('East [m]')
zlabel('Up [m]')
title('Interceptor and Target Trajectories')
legend('Interceptor','Target')
grid minor

%% Unreachable Set

%%% Target Trajectory

v2           = -2401;       % [m/s] Target Velocity - Mach 7
theta2       = 5;           % [deg] Flight path angle

t           = (0:.001:133)'; % time [s]


y2           = v2*sind(theta2)*t + h0;
x2           = v2*cosd(theta2)*t + x0;

figure
plot(t,sqrt((x2 - xI).^2 + (y2 - yI).^2))
xlabel('Time [s]')
ylabel('Distance [m]')
title('Difference in Position')
grid minor

figure
plot3(t,xI,yI,t,x2,y2)
zlim([0 max(y2)]);
xlabel('Time [s]')
ylabel('East [m]')
zlabel('Up [m]')
title('Interceptor and Target Trajectories')
legend('Interceptor','Target')
grid minor
