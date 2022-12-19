close all
clear
clc

%% Setup Measurements

% Truth Stae
xtrue       = [1 1 1]';

% Time Vector
time        = (0:.1:1)';

% Truth
y           = xtrue(1) + xtrue(2)*sin(10.*time) + xtrue(3)*exp(2*time.^2);

% Measurement error standard deviation
sigma       = [.001, .002, .005, .01, .008, .002,.01,.007,.02,.006,.001]';

% Measurement Error Covariance Matrix
R           = diag(sigma.^2);       % E{vvT} = R

% Measurments
ytilde      = y + sqrt(diag(R)).*(randn(length(time),1));

%% Minimum Variance with A Priori Estimate

% a priori x estimates
xa          = [1.01, 0.98, 0.99]';

% a priori covarince matrix
Q           = 0.001*eye(3);

% Basis Functions
H           = [ones(length(time),1), sin(10.*time), exp(2*time.^2)];

% Estimate Covarince Matrix
P           = inv(H'*inv(R)*H + inv(Q))

% Minimum Variance Estimate
xhat        = P*(H'*inv(R)*ytilde + inv(Q)*xa)

%% Residuals
r           = ytilde - (xhat(1) + xhat(2)*sin(10.*time) + xhat(3)*exp(2*time.^2));

% mean of residuals
resMean     = sum(r)/length(r)

% standard deviation of residuals
resStdDev   = sqrt(sum(r.^2)/(length(r) - 1))





