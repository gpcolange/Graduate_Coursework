clear
close all
clc

% Time
dt      = 1;
tend    = 3600;
t       = (0:dt:tend)';

% True Angular Rate
wtrue   = 0.0011;

% Gyro and Attitude Noise Standard Deviations
sigu    = sqrt(10)*1e-10;
sigv    = sqrt(10)*1e-7;
sign    = 17*1e-6;

% Attitude Measurements
ym      = t*wtrue + sign*randn(length(t),1);

% Gyro Transfer function
num_g   = dt*[1 1];
den_g   = 2*[1 -1];

% Gyro State Space Model
[phi_g,gam_g,c_g,d_g]   = tf2ss(num_g,den_g);

% Discrete Bias Propogation
bias    = dlsim(phi_g,gam_g,c_g,d_g,sigu/sqrt(dt)*randn(length(t),1),0.1*pi/180/3600/dt);

% Gyro Measurments
wm      = wtrue+sqrt(sigv^2/dt+1/12*sigu^2*dt)*randn(length(t),1)+bias;

% Discrete-Time Process Noise Covriance 
q       = [sigv^2*dt+1/3*sigu^2*dt^3 -1/2*sigu^2*dt^2;-1/2*sigu^2*dt^2 sigu^2*dt];

% Discrete State Space Matrices
phi     = [1 -dt;0 1];
gam     = [dt;0];

% Initial Covariance
poa         = 1e-4;
pog         = 1e-12;
p           = [poa 0;0 pog];
pcov        = zeros(length(t),2);
pcov(1,:)   = [poa pog];

% Initial Condition and H Matrix (constant)
x0          = [ym(1);0];
xe          = zeros(length(t),2);
xe(1,:)     = x0';
x           = x0;
h           = [1 0];

% Main Loop
for i = 1:length(t)-1

% Kalman Gain 
gain    = p*h'*inv(h*p*h'+sign^2);

% Update state and covariance
x       = x+gain*(ym(i)-h*x);
p       =[eye(2)-gain*h]*p;

% Propagate
x       = phi*x+gam*wm(i);
p       = phi*p*phi'+q;

% Store Variables
xe(i+1,:)   = x';
pcov(i+1,:) = diag(p)';

end

% 3-Sigma Outlier
sig3    = 3*sqrt(pcov);

% Plot Results
figure
% plot(t/60,[sig3(:,1) xe(:,1)-t*wtrue -sig3(:,1)]*1e6)
plot(t/60, sig3(:,1)*1e6, '--r', t/60, -sig3(:,1)*1e6, '--r', t/60,(xe(:,1)-t*wtrue)*1e6)
grid minor
xlabel('Time (Min)');
ylabel(' {Attitude Error ({\mu}rad)}');

figure
plot(t/60,xe(:,2)*180*3600/pi);grid
xlabel('Time (Min)');
ylabel('Bias Estimate (Deg/Hr)')