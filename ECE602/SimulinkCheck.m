HW7
out = sim('HW7_simulink.slx',time);

figure
subplot(4,1,1)
sgtitle('Part 3: Compare')
plot(T,X_lin(:,1),out.linear.time,out.linear.signals.values(:,1),'--k')
legend('MATLAB','Simulink')
ylabel('x [m]')
grid minor
subplot(4,1,2)
plot(T,X_lin(:,2),out.linear.time,out.linear.signals.values(:,2),'--k')
ylabel('$\dot{x}$ [m/s]','Interpreter','latex')
grid minor
subplot(4,1,3)
plot(T,X_lin(:,3),out.linear.time,out.linear.signals.values(:,3),'--k')
grid minor
ylabel('\theta [rad]')
subplot(4,1,4)
plot(T,X_lin(:,4),out.linear.time,out.linear.signals.values(:,4),'--k')
grid minor
ylabel('$\dot{\theta}$ [rad/s]','Interpreter','latex')
xlabel('Time [s]')

figure
subplot(4,1,1)
sgtitle('Part 4: Compare')
plot(out.nonlinear.time,out.nonlinear.signals.values(:,1) - X(:,1))
legend('Simulink','MATLAB')
ylabel('x [m]')
grid minor
subplot(4,1,2)
plot(out.nonlinear.time,out.nonlinear.signals.values(:,2) - X(:,2))
ylabel('$\dot{x}$ [m/s]','Interpreter','latex')
grid minor
subplot(4,1,3)
plot(out.nonlinear.time,out.nonlinear.signals.values(:,3) - X(:,3))
grid minor
ylabel('\theta [rad]')
subplot(4,1,4)
plot(out.nonlinear.time,out.nonlinear.signals.values(:,4) - X(:,4))
grid minor
ylabel('$\dot{\theta}$ [rad/s]','Interpreter','latex')
xlabel('Time [s]')

figure
subplot(4,1,1)
sgtitle('Part 5: Compare')
plot(T,X_lin_obs(:,1) - out.linearObs.signals.values(:,1),out.linearObs.time,out.linearObs.signals.values(:,3) - X_lin_obs(:,5))
legend('Difference between states','Difference between estimates')
ylabel('x [m]')
grid minor
subplot(4,1,2)
plot(T,X_lin_obs(:,2),T,X_lin_obs(:,6),out.linearObs.time,out.linearObs.signals.values(:,4),'--k')
ylabel('$\dot{x}$ [m/s]','Interpreter','latex')
grid minor
subplot(4,1,3)
plot(T,X_lin_obs(:,3) - out.linearObs.signals.values(:,2),out.linearObs.time,out.linearObs.signals.values(:,5) - X_lin_obs(:,7))
grid minor
ylabel('\theta [rad]')
subplot(4,1,4)
plot(T,X_lin_obs(:,4),T,X_lin_obs(:,8),out.linearObs.time,out.linearObs.signals.values(:,6),'--k')
grid minor
ylabel('$\dot{\theta}$ [rad/s]','Interpreter','latex')
xlabel('Time [s]')

figure
subplot(4,1,1)
sgtitle('Part 7: Compare')
plot(T,X_obs(:,1) - out.NonlinearObs.signals.values(:,1),out.NonlinearObs.time,out.NonlinearObs.signals.values(:,3) - X_obs(:,5))
legend('Difference between states','Difference between estimates')
ylabel('x [m]')
grid minor
subplot(4,1,2)
plot(T,X_obs(:,2),T,X_obs(:,6),out.NonlinearObs.time,out.NonlinearObs.signals.values(:,4),'--k')
ylabel('$\dot{x}$ [m/s]','Interpreter','latex')
grid minor
subplot(4,1,3)
plot(T,X_obs(:,3) - out.NonlinearObs.signals.values(:,2),out.NonlinearObs.time,out.NonlinearObs.signals.values(:,5) - X_obs(:,7))
grid minor
ylabel('\theta [rad]')
subplot(4,1,4)
plot(T,X_obs(:,4),T,X_obs(:,8),out.NonlinearObs.time,out.NonlinearObs.signals.values(:,6),'--k')
grid minor
ylabel('$\dot{\theta}$ [rad/s]','Interpreter','latex')
xlabel('Time [s]')