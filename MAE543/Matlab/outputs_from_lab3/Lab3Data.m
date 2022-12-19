clear
close all
clc

data        = importdata('sinusoid_0.8Hz_1V.csv');

t           = data.data(:,1)*.01;
theta_sim   = data.data(:,2);
theta_meas  = data.data(:,4);
input       = data.data(:,6);

figure
plot(t,theta_sim,t,theta_meas,t,input)
xlabel('Time [s]')
ylabel('Magnitude')
legend('Simulated Motor Position [rad]','Measured Motor Position [rad]','Input Voltage [V]')
grid minor

figure
plot(t,theta_sim - theta_meas)
xlabel('Time [s]')
ylabel('Error [rad]')
title('Position Error')
grid minor

