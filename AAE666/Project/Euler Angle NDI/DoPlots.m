figure
subplot(321)
plot(tout,phic*180/pi,tout,phic_filt*180/pi,tout,phi_p*180/pi,'--k')
grid minor
legend('Raw Command','Filtered Command','Response')
ylabel('Roll [deg]')
subplot(322)
plot(tout,phic_filt*180/pi - phi_p*180/pi)
grid minor
ylabel('Roll Error [deg]')
title('Residual')
subplot(323)
plot(tout,thetac*180/pi,tout,thetac_filt*180/pi,tout,theta_p*180/pi,'--k')
grid minor
ylabel('Pitch [deg]')
subplot(324)
plot(tout,thetac_filt*180/pi - theta_p*180/pi)
grid minor
ylabel('Pitch Error [deg]')
subplot(325)
plot(tout,psic*180/pi,tout,psic_filt*180/pi,tout,psi_p*180/pi,'--k')
grid minor
ylabel('Yaw [deg]')
xlabel('Time [s]')
subplot(326)
plot(tout,psic_filt*180/pi - psi_p*180/pi)
grid minor
ylabel('Yaw Error [rad]')
sgtitle(['Euler Angles: ',title_str])

figure
subplot(311)
plot(tout,w1*180/pi)
grid minor
ylabel('\omega_1 [deg/s]')
subplot(312)
plot(tout,w2*180/pi)
grid minor
ylabel('\omega_2 [deg/s]')
subplot(313)
plot(tout,w3*180/pi)
grid minor
ylabel('\omega_3 [deg/s]')
xlabel('Time [s]')
sgtitle(['Angular Velocity: ',title_str])

if dist ~= 0
    figure
    plot(tout,disturbance)
    xlabel('Time [s]')
    ylabel('Torque [N-m]')
    grid minor
    title('Disturbance torque applied to all 3 axis')
end

figure
plot(tout,u)
legend('u1','u2','u3')
title(['Control Torques: ', title_str])
ylabel('[N-m]')
xlabel('Time [s]')
grid minor