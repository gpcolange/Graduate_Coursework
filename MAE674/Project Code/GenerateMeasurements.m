GenerateTruth

%% Gyroscope Measurements
% Gyroscope Standard Deviation
sigu        = sqrt(10)*1e-10;   % rad/s^(3/2)
sigv        = sqrt(10)*1e-7;    % rad/s^(1/2)  

% Gyroscope Spectral Densities
Qu          = sigu^2 * eye(3);
Qv          = sigv^2 * eye(3);

% Initial Gyro Drift/Bias [rad/s]
b0          = 0.1*pi/180/3600*ones(3,1);

Gyro        = GyroscopeMeasurement(dt,Qu,Qv, b0, w_B_BI,t);

figure
subplot(3,1,1)
plot(t,Gyro.Measurements(1,:))
ylabel('Angular Rate [rad/s]')
title('Gyroscope Measurment $\hat{b}_1$','Interpreter','latex')
grid minor

subplot(3,1,2)
plot(t,Gyro.Measurements(2,:))
ylabel('Angular Rate [rad/s]')
title('Gyroscope Measurment $\hat{b}_2$','Interpreter','latex')
grid minor

subplot(3,1,3)
plot(t,Gyro.Measurements(3,:))
ylabel('Angular Rate [rad/s]')
title('Gyroscope Measurment $\hat{b}_3$','Interpreter','latex')
grid minor
xlabel('Time [s]')

figure
subplot(3,1,1)
plot(t,rad2deg(Gyro.Bias(1,:))*3600)
ylabel('Bias [deg/h]')
title('Axis 1 Gyro Drift')
grid minor
subplot(3,1,2)
plot(t,rad2deg(Gyro.Bias(2,:))*3600)
ylabel('Bias [deg/h]')
title('Axis 2 Gyro Drift')
grid minor
subplot(3,1,3)
plot(t,rad2deg(Gyro.Bias(3,:))*3600)
ylabel('Bias [deg/h]')
title('Axis 3 Gyro Drift')
grid minor
xlabel('Time [s]')

%% Star Tracker Measurments

% Angles Degrees
ascenion        = [0, 15, 30, 45, 60]';
declination     = [0, 30, 45, 60, 75]'; 

% Inertial Vectors
r1              = [cosd(declination(1))*cosd(ascenion(1));...
                   cosd(declination(1))*sind(ascenion(1));...
                   sind(declination(1))].*ones(3,length(t));
               
r2              = [cosd(declination(2))*cosd(ascenion(2));...
                   cosd(declination(2))*sind(ascenion(2));...
                   sind(declination(2))].*ones(3,length(t));
               
r3              = [cosd(declination(3))*cosd(ascenion(3));...
                   cosd(declination(3))*sind(ascenion(3));...
                   sind(declination(3))].*ones(3,length(t));
               
r4              = [cosd(declination(4))*cosd(ascenion(4));...
                   cosd(declination(4))*sind(ascenion(4));...
                   sind(declination(4))].*ones(3,length(t));
               
r5              = [cosd(declination(5))*cosd(ascenion(5));...
                   cosd(declination(5))*sind(ascenion(5));...
                   sind(declination(5))].*ones(3,length(t));

% Sensor Standard Deviations [rad]
sig1            =  3e-5;
sig2            =  3e-5;
sig3            =  3e-5;
sig4            =  3e-5;
sig5            =  3e-5;

R1              = sig1^2*eye(3);
R2              = sig2^2*eye(3);
R3              = sig3^2*eye(3);
R4              = sig4^2*eye(3);
R5              = sig5^2*eye(3);

StarTracker1    = StarTrackerMeasurements(q, R1, r1);
StarTracker2    = StarTrackerMeasurements(q, R2, r2); 
StarTracker3    = StarTrackerMeasurements(q, R3, r3);
StarTracker4    = StarTrackerMeasurements(q, R4, r4);
StarTracker5    = StarTrackerMeasurements(q, R5, r5);

% Measurement Error Covariance Matrix
R               = diag([diag(R1); diag(R2); diag(R3); diag(R4);diag(R5)]);

% Star Observations
btilde1         = StarTracker1.BodyMeasurements;
btilde2         = StarTracker2.BodyMeasurements;
btilde3         = StarTracker3.BodyMeasurements;
btilde4         = StarTracker4.BodyMeasurements;
btilde5         = StarTracker5.BodyMeasurements;

btilde          = [btilde1; btilde2; btilde3; btilde4;btilde5];

figure
subplot(3,1,1)
plot(t,btilde1(1,:),t,StarTracker1.BodyTrue(1,:))
grid minor
ylabel('$b_x \hat{b}_1$ [m] ', 'Interpreter','latex')
title('Star Observations')
legend('$\tilde{b}_1$','$b_1$', 'Interpreter','latex')
subplot(3,1,2)
plot(t,btilde1(2,:),t,StarTracker1.BodyTrue(2,:))
grid minor
ylabel('$b_y \hat{b}_2$ [m] ', 'Interpreter','latex')
subplot(3,1,3)
plot(t,btilde1(3,:),t,StarTracker1.BodyTrue(3,:))
grid minor
ylabel('$b_z \hat{b}_3$ [m] ', 'Interpreter','latex')
xlabel('Time [s]')

