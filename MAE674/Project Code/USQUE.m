%% Unscented Quaternion Estimator (USQUE)

%% Define and Initialize 
% Measurment Error Covariance
R               = diag([diag(sig1^2*eye(3));diag(sig2^2*eye(3));diag(sig3^2*eye(3));...
                        diag(sig4^2*eye(3));diag(sig5^2*eye(3))]);

% Local Error Quaternion scale factor
a               = 1;
f               = 2*(a+1);

% Initial State Error Covariance Matrix
P0              = diag([P0a;P0b]);
P               = P0;
UKF.Pdiag       = zeros(length(P0),length(t));
UKF.Pdiag(:,1)  = diag(P);

% Initial State Estimates
delp0           = zeros(3,1);
UKF.xhat        = zeros(length(delp0)+length(Bhat0),length(t));
UKF.xhat(:,1)   = [delp0;Bhat0];

% Initial Quaternion Estimate
UKF.qhat        = zeros(length(qhat0),length(t));
UKF.qhat(:,1)   = qhat0;

% Calculate Qbar - Equation 43
Qbar            = (dt/2)*[(sigv^2 - sigu^2*(1/6)*dt^2)*eye(3), zeros(3);...
                     zeros(3), sigu^2*eye(3)];

% UKF Parameters
lambda          = 5;  % Chosen tuning parameter
n               = length(UKF.xhat(:,1));

%% Run UKF From Paper
for i = 1:length(t)-1

    % Equation 5 - Calculate Sigma Points
    squareroot  = chol((n+lambda)*(P + Qbar))';
    sigk        = real([squareroot, -squareroot]);
    % X_k(0)
    Chik0       = UKF.xhat(:,i);    
    % X_k(i) -> i = 1:12
    Chiki       = sigk + kron(UKF.xhat(:,i),ones(1,2*n));

    % Partition Chi0 - Equation 31
    % X_k(0) = UKF.xhat_k+ = [del_pk+, Beta_k+]'
    delpk       = Chik0(1:3);
    Betahat     = Chik0(3:end);

    % Equation 32
    % X_k(i) = [X_delp_k(i);X_Beta_k(i)]' -> i = 0:12
    Chik        = [Chik0,Chiki];
    Chi_delpi   = Chik(1:3,:);
    Chi_betai   = Chik(4:6,:);

    % Equation 33a
    % UKF.qhat_k(0)
    qhatk0      = UKF.qhat(:,i);
    % UKF.qhat_k(i) -> i = 1:12
    for j = 1:12

        % Equation 34 - i = 1:12... ignore j = 1 thats my 0 term
        delq4i      = -a*norm(Chi_delpi(:,j+1))^2 + f*sqrt(f^2 + (1 - a^2)*norm(Chi_delpi(:,j+1))^2)/...
                        (f^2 + norm(Chi_delpi(:,j+1))^2);
        delrhoi     = (f^-1)*(a + delq4i)*Chi_delpi(:,j+1);

        % delqk+(i) = [delrho+'(i), delq4+(i)]'
        delqki(:,j) = [delrhoi' delq4i]';

        % Equation 33b UKF.qhatk+(i) - i =  1:12
        qhatki(:,j) = QuaternionMultiply(delqki(:,j),UKF.qhat(:,i));
    end

    % Propogate Quaternion 
    % qhat+_k
    qhatk   = [qhatk0,qhatki];

    for j = 1:13
        % Equation 36 i = 0:12
        omegahat   = Gyro.Measurements(:,i) - Chi_betai(:,j);

        qhatkp1i(:,j)   = RungeKutta4(@QuaternionKinematics,t(i),qhatk(:,j),dt,omegahat);

        % [UKF.qhat-_k+1(0)]^-1
        invqhatkp10     = [-qhatkp1i(1,1), -qhatkp1i(2,1), -qhatkp1i(3,1), qhatkp1i(4,1)]';

        % Equation 37 delq-_k+1(i) i = 0:12
        delqkp1(:,j)    =  QuaternionMultiply(qhatkp1i(:,j),invqhatkp10);
    end

    % Propogate Sigma Points: X_delp_k+1(0), X_delp_k+1(i)
    Chi_delp_kp10       = zeros(3,1);

    for j  = 1:12
        % Equation 38
        delrhokp1           = delqkp1(1:3,j+1);
        delq4kp1            = delqkp1(4,j+1);

        Chi_delp_kp1i(:,j)  = f*delrhokp1/(a + delq4kp1);
    end

    %X_K+1(0)
    Chi_kp10    = [Chi_delp_kp10;Chi_betai(:,1)];

    %X_k+1(i)
    Chi_kp1i    = [Chi_delp_kp1i;Chi_betai(:,2:end)];

    % Predicted mean equation 7 - UKF.xhat-_k+1
    xhatkp1     = (1/(n+lambda))*(lambda*Chi_kp10 + (1/2)*sum(Chi_kp1i,2));

    % Predicted Covariance equation 8 - P-_k+1
    Pmat        = Chi_kp1i - kron(xhatkp1,ones(1,2*n));
    Pkp1        = (1/(n+lambda))*(lambda*(Chi_kp10 - xhatkp1)*(Chi_kp10 - xhatkp1)'...
                   + (1/2)*Pmat*Pmat')+ Qbar;

    % Observations Equation 44
    for j  = 1:13
        Ahat        = Quaternion2DCM(qhatkp1i(:,j));
        gamma(:,j)  = [Ahat*r1(:,i+1);Ahat*r2(:,i+1);Ahat*r3(:,i+1);Ahat*r4(:,i+1);Ahat*r5(:,i+1)];
    end

    % mean observation equation 9
    yhat        = (1/(n+lambda))*(lambda*gamma(:,1) + (1/2)*sum(gamma(:,2:end),2));

    % Output Covariance Pyy_k+1 Equation 11
    Pyysum      = lambda*(gamma(:,1) - yhat)*(gamma(:,1) - yhat)';

    for j = 1:(2*n)
        Pyysum = Pyysum + (1/2)*(gamma(:,j+1) - yhat)*(gamma(:,j+1) - yhat)';
    end

    Pyy         = Pyysum*(1/(n+lambda));

    % Innovation Covariance Equation 12
    Pvv         = Pyy + R;

    % Cross Correlation Equation 13
    Pxysum      = lambda*(Chi_kp10 - xhatkp1)*(gamma(:,1) - yhat)';

    for j = 1:(2*n)
        Pxysum = Pxysum + (1/2)*(Chi_kp1i(:,j) - xhatkp1)*(gamma(:,j+1) - yhat)';
    end

    Pxy         = Pxysum*(1/(n+lambda));
  
    % innovation equation 3
    vk          = btilde(:,i+1) - yhat;

    % Kalman Gain equation 4
    K           = real(Pxy*inv(Pvv));

    % State and Covariance Update equation 2
    UKF.xhat(:,i+1) = xhatkp1 + K*vk;

    P               = Pkp1 - K*Pvv*K';
    UKF.Pdiag(:,i+1)= diag(P);

    % Equation 45
    delp         = UKF.xhat(1:3,i+1);
    delq4        = -a*norm(delp)^2 + f*sqrt(f^2 + (1 - a^2)*norm(delp)^2)...
                    /(f^2 + norm(delp)^2);
    delrho       = f^-1 * (a + delq4)*delp;
    delq         = [delrho;delq4];

    % Equation 44
    UKF.qhat(:,i+1)  = QuaternionMultiply(delq,qhatkp1i(:,1));
    
    if abs(1 - norm(UKF.qhat(:,i+1) )) > 1e-7
        UKF.qhat(:,i+1) = UKF.qhat(:,i+1)/norm(UKF.qhat(:,i+1));
    end

    delp         = zeros(size(delp));

end

%% Calculate Small Angle Errors

% Initialize Error quaternion
delqtrue        = zeros(length(qhat0), length(t));

% Initialize Small Yaw Pitch Roll Errors
UKF.roll_error  = zeros(1,length(t));
UKF.pitch_error = zeros(1,length(t));
UKF.yaw_error   = zeros(1,length(t));

for i = 1:length(t)
    % Compute Error Quaternion: delq = q x UKF.qhat^-1 
    % function calculates UKF.qhat^-1 from UKF.qhat
    delqtrue(:,i)       = QuaternionError(q(:,i),UKF.qhat(:,i));
    
    % Calculate Yaw, Pitch, and Roll Attitude Errors
    att_err             = delqtrue(1:3,i)*2;
    UKF.roll_error(i)   = rad2deg(att_err(1));
    UKF.pitch_error(i)  = rad2deg(att_err(2));
    UKF.yaw_error(i)    = rad2deg(att_err(3));

end



    