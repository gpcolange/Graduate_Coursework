clear
close all
clc

data        = importdata('PD_square_0.4Hz_0.5_Kp10_Kd_0.1.csv');
vdata       = importdata('VoltagePD_square_0.4Hz_0.5_Kp10_Kd_0.1.csv');

t           = data.data(:,1)*.01;
theta_sim   = data.data(:,2);
theta_meas  = data.data(:,4);
input       = data.data(:,6);
tv          = vdata.data(:,1)*.01;
vsim        = vdata.data(:,2);
v           = vdata.data(:,4);

figure
plot(t,theta_sim,t,theta_meas,t,input)
xlabel('Time [s]')
ylabel('Motor Position [rad]')
legend('Simulated','Measured','Reference')
title('Motor Positions')
grid minor

figure
plot(tv,vsim, t, v)
xlabel('Time [s]')
ylabel('Voltage[V]')
legend('Simulated','Actual')
title('Control Input Voltage')
grid minor
