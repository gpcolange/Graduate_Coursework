%% MAE 562 Research Paper - Motion Camouflague

clear
close all
clc


mu          = 10                                                        ;   % Control Gain
r0          = [1500, 500]'                                              ;   % Initial Distance [m]
rdot0       = [0 0]'                                                    ;   % Intial velocity [m/s]

time        = (0:.05:10)'                                               ;   % time [s]
IC          = [r0 rdot0]                                                ;   % Initial Condition Vector

options     = odeset('AbsTol',1e-3,'RelTol',1e-3,'Events', @BREAK)      ;   % ODE45 solver options
[T,Z]       = ode45(@(t,z)MC(t,z,mu),time,IC,options)                   ;

x           = Z(:,1);
y           = Z(:,2);
xdot        = Z(:,3);
ydot        = Z(:,4);

r           = sqrt(x.^2 + y.^2);
rdot        = sqrt(xdot.^2 + ydot.^2);

figure
ax1         = subplot(2,1,1);
plot(T,r)
xlabel('Time [s]')
ylabel('$\left|r\right|$ [m]','Interpreter','Latex')
grid minor
title('Motion Camouflage')

ax2         = subplot(2,1,2);
plot(T,rdot)
xlabel('Time [s]')
ylabel(' $\left|\dot{r}\right|$ [m/s]','Interpreter','Latex')
grid minor

