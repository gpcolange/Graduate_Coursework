%% Sequential Least Squares vs Batch Least Squares

clear
clc
close all

%% Generate Truth Measurements from Model
% Discrete Model: y_k+1 = Phi*y_k + Gamma*u_k
% Continous Model: xdot = Ax + Bu
% Let A  = -1 and B = 1;

T               = 0.1;                                                      % Period [s]
time            = (0:T:10)';                                                % Time Vector
    
A               = -1;                                                       % Contious A matrix
B               = 1;                                                        % Continous B matrix

Phi             = expm(A*T);                                                % Discrete A matrix
Gamma           = (Phi - eye(size(Phi)))*inv(A)*B;                          % Discrete B matrix
Cd              = eye(size(Phi));                                           % Discrete C, full state feedback
Dd              = zeros(size(Gamma));                                       % Discrete D

u               = [100; zeros(length(time)-1,1)];                           % Impulse input
y               = dlsim(Phi,Gamma,Cd,Dd,u);                                 % Discrete Sim - Truth

std             = 0.08;                                                     % Measurement Noise standard deviation
ym              = y + std*randn(length(y),1);                               % Measurment

%% Batch Least Squares 
H               = [ym(1:end-1), u(1:end-1)];                                % H matrix, basis functions
xhat            = inv(H'*H)*H'*ym(2:end);                                   % Least squares estimate
PhiHat          = xhat(1);                                                  % Estiamted Phi
GammaHat        = xhat(2);                                                  % Estimated Gamma

yhat            = dlsim(xhat(1),xhat(2),Cd,Dd,u);                           % Estimated output

figure
plot(time,yhat,time,ym,'*')
xlabel('Time [s]')
ylabel('Output')
legend('Estimate','Measurment')
title('Standard Batch Least Squares')
grid minor

%% Sequential Least Squares

W               = std^-2;                                                   % Optimal Weight Matrix
alpha           = 1e3;                                                      % Initialization Parameter
beta            = [1e-2, 1e-2]';                                            % Initialization Parameter

% Initialize P and xhat
Pdiag           = zeros(2,length(ym)-1);
xhat_sls        = zeros(2,length(ym)-1);
yhat_sls        = zeros(length(ym)-1,1);

P1              = inv(1/(alpha^2)*eye(2) + H(1,:)'*W*H(1,:));
xhat_sls(:,1)   = (P1*(beta/alpha + H(1,:)'*W*ym(2,:)));

Pdiag(:,1)      = diag(P1);

Pk              = P1;


for i = 1:length(ym)-2
    K               = Pk*H(i+1,:)'*inv(H(i+1,:)*Pk*H(i+1,:)' + inv(W)); % K_k+1
    xhat_sls(:,i+1) = xhat_sls(:,i) + K*(ym(i+2,:) - H(i+1,:)*xhat_sls(:,i)); % x_k+1
    Pk1             = (eye(size(P1)) - K*H(i+1,:))*Pk; % P_k+1
    Pdiag(:,i+1)    = diag(Pk1);
    Pk              = Pk1;
end

figure
ax1 = subplot(2,1,1);
plot(time(1:end-1),xhat_sls(1,:))
ylim([0 10]);
ylabel('Magnitude - \Phi')
legend('\Phi')
grid minor

ax2 = subplot(2,1,2);
plot(time(1:end-1),xhat_sls(2,:))
ylim([xhat_sls(2,end)-.001, xhat_sls(2,end)+.001])
ylabel('Magnitude - \Gamma')
legend('\Gamma')
xlabel('Time [s]')
grid minor
linkaxes([ax1, ax2],'x')
sgtitle('xhat -  Sequential Least Squares')


figure
h1 = subplot(2,1,1);
plot(time(1:end-1),Pdiag(1,:))
ylabel('Magnitude')
legend('P11')
grid minor

h2 = subplot(2,1,2);
plot(time(1:end-1),Pdiag(2,:))
ylabel('Magnitude')
legend('P22')
grid minor
xlabel('Time [s]')
linkaxes([h1, h2],'x')
sgtitle('P diagonal')

yhat_sls(1) = 0;

% Calculate yhat using current estimate of phi and gamma 
for i = 1:length(xhat_sls)-1
    yhat_sls(i+1) = xhat_sls(1,i)*yhat_sls(i) + xhat_sls(2,i)*u(i);
end

figure
plot(time(1:end-1),yhat_sls,time,ym,'*',time,yhat)
ylabel('Output')
legend('Sequential Estimate','Measurment','Batch Estimate')
grid minor
title('Sequential vs Least Squares')

%% Crassidis Code for Sequential Least Squares
% Weight and H Matrix
dt=0.1;tf=10;
t=[0:dt:tf];
m=length(t);
w=inv(0.08^2);
h=[ym(1:m-1) u(1:m-1)];

% Initial Conditions for Sequential Algorithm
alpha=1e3;
beta=[1e-2;1e-2];
p0=inv(1/alpha/alpha*eye(2)+h(1,:)'*w*h(1,:));
x0=p0*(1/alpha*beta+h(1,:)'*w*ym(2));

% Sequential Least Squares
xr=zeros(m-1,2);xr(1,:)=x0';
p=zeros(m-1,2);p(1,:)=diag(p0)';pp=p0;
for i=1:m-2;
 k=pp*h(i+1,:)'*inv(h(i+1,:)*pp*h(i+1,:)'+inv(w));
 pp=(eye(2)-k*h(i+1,:))*pp;
 xr(i+1,:)=xr(i,:)+(k*(ym(i+2)-h(i+1,:)*xr(i,:)'))';
 p(i+1,:)=diag(pp)';
end

% Plot Results
figure
subplot(221)
plot(t(1:100),xr(:,1));grid
axis([0 10 0 10]);
set(gca,'fontsize',12);
set(gca,'xtick',[0 2 4 6 8 10]);
set(gca,'ytick',[0 2 4 6 8 10]);
xlabel('Time (Sec)')
ylabel('{\it x}_1 Estimate')

subplot(222)
plot(t(1:100),xr(:,2));grid
axis([0 10 0.0945 0.096]);
set(gca,'fontsize',12);
set(gca,'xtick',[0 2 4 6 8 10]);
set(gca,'ytick',[0.0945 0.095 0.0955 0.096]);
xlabel('Time (Sec)')
ylabel('{\it x}_2 Estimate')

subplot(223)
semilogy(t(1:100),p(:,1));grid
axis([0 10 1e-5 1e5])
set(gca,'xtick',[0 2 4 6 8 10]);
set(gca,'ytick',[1e-5 1 1e5])
xlabel('Time (Sec)')
ylabel('{\it P}_{11}')

subplot(224)
semilogy(t(1:100),p(:,2));grid
set(gca,'xtick',[0 2 4 6 8 10]);
xlabel('Time (Sec)')
ylabel('{\it P}_{22}')
sgtitle('Crassidis Sequential Least Squares Results')
