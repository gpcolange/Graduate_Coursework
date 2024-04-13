figure
subplot(411)
plot(tout,q1c(:),tout,q1,'--r')
grid minor
legend('Command','Response')
ylabel('q_1')
subplot(412)
plot(tout,q2c(:),tout,q2,'--r')
grid minor
ylabel('q_2')
subplot(413)
plot(tout,q3c(:),tout,q3,'--r')
grid minor
ylabel('q_3')
xlabel('Time [s]')
sgtitle(title_str1)
subplot(414)
plot(tout,vecnorm([q1c(:)';q2c(:)';q3c(:)';q4c(:)']),tout,vecnorm([q1';q2';q3';q4']),'--r')
ylabel('||q||')
title("Quaternion Norm")
grid minor
xlabel('Time [s]')

% Convert to Euler Angles from Quaternions
[phi_calc,theta_calc, psi_calc] = QuaterniontoEuler321(q1,q2,q3,q4);

figure
subplot(311)
plot(tout,phic,tout,phic_filt,tout,phi_calc,'--k')
grid minor
ylabel('Roll [rad]')
legend('Raw Command','Filtered Command','Response')
subplot(312)
plot(tout,thetac,tout,thetac_filt,tout,theta_calc,'--k')
grid minor
ylabel('Pitch [rad]')
subplot(313)
plot(tout,psic,tout,psic_filt,tout,psi_calc,'--k')
grid minor
ylabel('Yaw [rad]')
xlabel('Time [s]')
sgtitle(title_str2)
