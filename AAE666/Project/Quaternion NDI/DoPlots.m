
for i = 1:length(tout)
    dq(:,i) = QuaternionError([q1c(i);q2c(i);q3c(i);q4c(i)],[q1(i);q2(i);q3(i);q4(i)]);
end

figure
subplot(421)
plot(tout,q1c(:),tout,q1(:),'--r')
grid minor
legend('Command','Response')
ylabel('q_1')
title('Response')
subplot(422)
plot(tout,dq(1,:))
grid minor
ylabel('\delta q_1')
title('Multiplicative Error')
subplot(423)
plot(tout,q2c(:),tout,q2(:),'--r')
grid minor
ylabel('q_2')
subplot(424)
plot(tout,dq(2,:))
grid minor
ylabel('\delta q_2')
subplot(425)
plot(tout,q3c(:),tout,q3(:),'--r')
grid minor
ylabel('q_3')
xlabel('Time [s]')
subplot(426)
plot(tout,dq(3,:))
grid minor
ylabel('\delta q_3')
sgtitle(['Quaternions: ',title_str])
subplot(427)
plot(tout,vecnorm([q1c(:)';q2c(:)';q3c(:)';q4c(:)']),tout,vecnorm([q1(:)';q2(:)';q3(:)';q4(:)']),'--r')
ylabel('||q||')
title("Quaternion Norm")
grid minor
xlabel('Time [s]')
subplot(428)
plot(tout,dq(4,:))
grid minor
ylabel('\delta q_4')

% Convert to Euler Angles from Quaternions
[phi_calc,theta_calc, psi_calc] = QuaterniontoEuler321(q1(:),q2(:),q3(:),q4(:));

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
sgtitle(['Euler Angles: ',title_str])

figure
subplot(311)
plot(tout,w1)
grid minor
ylabel('\omega_1 [rad/s]')
subplot(312)
plot(tout,w2)
grid minor
ylabel('\omega_2 [rad/s]')
subplot(313)
plot(tout,w3)
grid minor
ylabel('\omega_3 [rad/s]')
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
