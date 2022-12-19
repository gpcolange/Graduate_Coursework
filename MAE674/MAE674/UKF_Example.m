%% UKF Example Lecture 15

clear
close all
clc

%% Time
dt      = .1;
t       = 0:dt:60;

%% Truth Model

F       = [-4 -3 -4 -1;...
            1 0 0 0;...
            0 1 0 0;
            0 0 1 0];

H       = [1 0 0 0; 0 0 1 0];

%% Initial Conditions
x0      = [1 0 2 0]';

% Initial State Estimates
xhat0   = zeros(4,1);

% Initial Error Covariance
P0      = diag([(2/3)^2 .001 (4/3)^2 .001]);

[n,col] = size(F);

% Measurment Error Covariance
R       = [.01 .005; .005 .02];

% Continous Time Process Noise
Q       = .01;

G       = [1 0 0 0]';

% Convert to discrete
A       = [-F, G*Q*G';zeros(size(F)), F']*dt;
B       = expm(A);
Phi     = (B(n+1:end,n+1:end))';
Qk      = Phi*B(1:n,n+1:end);


%% Generate Correlated Noise

% Eigenvalue decompisition of R = VDV'
[V,D]    = eig(R);

v_uncor  = (randn(length(t),2)*sqrt(D))';
v        = V*v_uncor;

% Eigenvalue decomposition of Qk
[V,D]    = eig(Qk);
w_uncor  = (randn(length(t),n)*sqrt(D))';
w        = V*w_uncor;

%% Measurments
x           = zeros(4,length(t));
x(:,1)      = x0;

for i = 1:length(t)-1
    x(:,i+1) = Phi*x(:,i) + w(:,i);
end

y           = H*x;
ytilde      = y + [v(1,:);v(2,:)];

%% Unscented Kalman Filter

xhat        = zeros(size(x));
xhat(:,1)   = xhat0;
Pdiag       = zeros(n,length(t));
Pdiag(:,1)  = diag(P0);
P           = P0;

% Augmented Covariance Matrix P_a = [P Pxw; Pxw' Qk] if R appears linearly
% Pxw is typically zero. L is number of columns from augmented P_a
P_a         = blkdiag(P0,Qk);
[row,L]     = size(P_a);

alpha       = 1;
beta        = 0;
kap         = 3 + L;
lambda      = alpha^2*(L+kap) - L;
gamma       = sqrt(L + lambda);

W0mean      = lambda/(L+lambda);
W0cov       = lambda/(L + lambda) + (1 -alpha^2 + beta);
Wimean      = 1/(2*(L+lambda));

% Break up matrix square root - Qk constant, take out for loop
squarerootq = chol(Qk)';
sigw        = real([zeros(n,2*n), gamma*squarerootq, - gamma*squarerootq]);

for i = 1:length(t)-1
    squarerootp = chol(P)';
    sigv        = real([gamma*squarerootp, - gamma*squarerootp]);
    Xa0         = xhat(:,i);
    Xa          =[sigv+kron(xhat(:,i),ones(1,2*n)) repmat(Xa0,1,2*n)];

    % Propogate States and Covariance
    Xa0         = Phi*Xa0;
    Xa          = Phi*Xa + sigw;
    xhat(:,i+1) = W0mean*Xa0 + Wimean*sum(Xa,2);

    P0          = W0cov*(Xa0 - xhat(:,i+1))*(Xa0 - xhat(:,i+1))';
    Pk          = Xa - kron(xhat(:,i+1),ones(1,4*n));
    P           = P0 + Wimean*Pk*Pk';

    % Output
    yhatk       = H*Xa;
    yhat0       = H*Xa0;
    yhat        = W0mean*yhat0 + Wimean*sum(yhatk,2);

    % Pyy and Pxy
    Pyy0        = W0cov*(yhat0 - yhat)*(yhat0 - yhat)';
    Pyyk        = yhatk - kron(yhat,ones(1,4*n));
    Pyy         = Pyy0 + Wimean*Pyyk*Pyyk';

    Pxy0        = W0cov*(Xa0 - xhat(:,i+1))*(yhat0 - yhat)';
    Pxy         = Pxy0 + Wimean*Pk*Pyyk';

    % Innovations Covariance
    Pvv        = Pyy + R;

    % Gain
    K           = real(Pxy*inv(Pvv));
    xhat(:,i+1) = xhat(:,i+1) + K*(ytilde(:,i+1) - yhat);
    P           = P - K*Pvv*K';
    Pdiag(:,i+1)= diag(P);

end


figure
subplot(4,1,1)
plot(t,x(1,:) - xhat(1,:),t,3*sqrt(Pdiag(1,:)),'--r',t,-3*sqrt(Pdiag(1,:)),'--r')
xlabel('Time [s]')
grid minor

subplot(4,1,2)
plot(t,x(2,:) - xhat(2,:),t,3*sqrt(Pdiag(2,:)),'--r',t,-3*sqrt(Pdiag(2,:)),'--r')
xlabel('Time [s]')
grid minor

subplot(4,1,3)
plot(t,x(3,:) - xhat(3,:),t,3*sqrt(Pdiag(3,:)),'--r',t,-3*sqrt(Pdiag(3,:)),'--r')
xlabel('Time [s]')
grid minor

subplot(4,1,4)
plot(t,x(4,:) - xhat(4,:),t,3*sqrt(Pdiag(4,:)),'--r',t,-3*sqrt(Pdiag(4,:)),'--r')
xlabel('Time [s]')
grid minor