%% Multiplicative Extended Kalman Filter

%% Initialize 

% Initial angular rate estimate
omegahat    = zeros(size(w_B_BI));

% State Estimates - [q1 q2 q3 q4 B1 B2 B3]'
EKF.xhat        = zeros(length(qhat0)+length(Bhat0),length(t));
EKF.xhat(:,1)   = [qhat0;Bhat0];

% Initial State Error Covariance Matrix
P0          = diag([P0a;P0b]);
P           = P0;

% Initialize diagonal of error covariance matrix
EKF.Pdiag       = zeros(length(P0),length(t));
EKF.Pdiag(:,1)  = diag(P0);

% Discrete time version of G Matrix
Gamma       = [-eye(3), zeros(3); zeros(3), eye(3)];

% Discrete Process Noise Covariance
Qk          = [(sigv^2*dt + (1/3)*sigu^2*dt^3)*eye(3), ((1/2)*sigu^2*dt^2)*eye(3);...
               ((1/2)*sigu^2*dt^2)*eye(3), sigu^2*dt*eye(3)];

% Continous Process Noise Covariance 
Q           = [sigv^2*eye(3), zeros(3); zeros(3), sigu^2*eye(3)];

% Measurement Error Covaraince
R               = diag([diag(sig1^2*eye(3));diag(sig2^2*eye(3));diag(sig3^2*eye(3));...
                        diag(sig4^2*eye(3));diag(sig5^2*eye(3))]);

for i = 1:length(t)-1
    Bhat    = EKF.xhat(5:end,i);

    % Attitude Matrix from Quaternion Estimate
    qhat    = EKF.xhat(1:4,i);
    Ahat    = Quaternion2DCM(qhat);

    %% Gain

    % Sensitivity Matrix
    H      = [ skew(Ahat*r1(:,i)), zeros(3);...
               skew(Ahat*r2(:,i)), zeros(3);...
               skew(Ahat*r3(:,i)), zeros(3);...
               skew(Ahat*r4(:,i)), zeros(3);...
               skew(Ahat*r5(:,i)), zeros(3)];

    % Kalman Gain
    K       = P*H'*inv((H*P*H' + R));

    %% Update

    % Update State Error Covariance Matrix
    P               = (eye(size(P)) - K*H)*P;

    % estimated output
    h               = [Ahat*r1(:,i);Ahat*r2(:,i);Ahat*r3(:,i);Ahat*r4(:,i);Ahat*r5(:,i)];

    % Update state error - btilde is star tracker body measurements
    % deltaEKF.xhat = [del_alpha deltaB] 
    % where alpha is quaternion small angle approximation del_q = del_alpha
    deltaEKF.xhat       = K*(btilde(:,i) - h);    

    % Components of state error
    del_alpha       = deltaEKF.xhat(1:3);
    delta_B         = deltaEKF.xhat(4:end);

    % Vector part of quaternion
    rho             = qhat(1:3);
    % E(qhat)
    E               = [(qhat(4)*eye(3)+skew(rho));-rho'];

    % q+_k - quaternion estimate update at current time 
    qhat            = qhat + (1/2)*E*del_alpha;
    
    if abs(1 - norm(qhat)) > 1e-7
        qhat        = qhat/norm(qhat);
    end

    % B+_k - bias estimate update at current time
    Bhat            = Bhat + delta_B;

    % xe+_k - full state estime update at current time
    EKF.xhat(:,i)   = [qhat;Bhat];

    %% Propogate
    
    % w+_k - angular velocity estiame at current time
    omegahat(:,i)       = Gyro.Measurements(:,i) - Bhat;
    
    % qdot = 1/2*E(q+_k)*what_k : q-_k+1 -> quaternion estimate for next
    % time step
    EKF.xhat(1:4,i+1)   = RungeKutta4(@QuaternionKinematics,t(i),qhat,dt,omegahat(:,i));

    % B-_k+1 = B+_k -> bias estimate for next time step
    EKF.xhat(5:end,i+1) = Bhat;

    % Error Model Matrices: delta_xdot = F*delta_x + G*w
    F               = [-skew(omegahat(:,i)), -eye(3); zeros(3), zeros(3)];

    % Convert F from continous to discrete
    Phi             = expm(F*dt);

    % Use discrete time solution for error covariance
    P               = Phi*P*Phi' + Gamma*Qk*Gamma';
    EKF.Pdiag(:,i+1)    = diag(P);

    % P contains covariances for atttitude yaw, pitch,roll and bias
    % so has 6 rows instead of 7 rows as seen in full state vector

end

%% Calculate Small Angle Errors

% Initialize Error quaternion
delq        = zeros(length(qhat0), length(t));

% Initialize Small Yaw Pitch Roll Errors
EKF.roll_error  = zeros(1,length(t));
EKF.pitch_error = zeros(1,length(t));
EKF.yaw_error   = zeros(1,length(t));

for i = 1:length(t)
    % Compute Error Quaternion: delq = q x qhat^-1 
    % function calculates qhat^-1 from qhat
    qhat                = EKF.xhat(1:4,i);
    Bhat                = EKF.xhat(5:end,i);
    EKF.delq(:,i)       = QuaternionError(q(:,i),qhat);
    
%     qhatinv         = [-qhat(1),-qhat(2),-qhat(3),qhat(4)]';
%     dqtest(:,i)     = QuaternionMultiply(q(:,i),qhatinv);
    
    % Calculate Yaw, Pitch, and Roll Attitude Errors
    att_err             = EKF.delq(1:3,i)*2;
    EKF.roll_error(i)   = rad2deg(att_err(1));
    EKF.pitch_error(i)  = rad2deg(att_err(2));
    EKF.yaw_error(i)    = rad2deg(att_err(3));

end

 