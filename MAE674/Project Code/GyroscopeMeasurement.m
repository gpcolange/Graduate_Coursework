function [Gyro] = GyroscopeMeasurement(dt,Qu, Qv, b0, w_B_BI,t)
%%% Gyroscope Measurement Model 
%%% Inputs are time vector, sampling period, spectral densities matrices of gaussian 
%%% noise processes Nu and Nv, and initial bias/drift in rad/s, and the true angular velocity 
%%% Output struct containing gyro bias (or drift) and measurements for each axis
%%% Variables:
%%% dt      = sampling period
%%% t       = time vector
%%% Qu      = spectral density of Nu - sigu^2*I_3x3
%%% Qv      = spectral density of Nv - sigv^2*I_3x3
%%% b0      = initial bias/drift
%%% w_B_BI  = true angular rate of body wrt to inertial frame, in body
%%% coordinates. 3xm matrix

sigu                = sqrt(diag(Qu));
sigv                = sqrt(diag(Qv));

% Initialize size
bias                = zeros(3,length(w_B_BI));
wtilde_B_BI         = bias;

% Zero mean Gaussian white-noise process Nu and Nv
Nu                   = randn(3,length(w_B_BI));
Nv                   = randn(3,length(w_B_BI));

bias(:,1)            = b0;
wtilde_B_BI(:,1)     = w_B_BI(:,1) + bias(:,1) + sqrt(((sigv.^2)/(dt)) + (1/12)*(sigu.^2)*dt).*Nv(:,1);

for j = 1:length(t)-1
    % b_k+1 = b_k + sigu*sqrt(dt)*Nu 
    bias(:,j+1)           = bias(:,j) + sigu*sqrt(dt).*Nu(:,j);
    
    % w_measured_k+1 = w_k+1 + 1/2*(b_k+1 + b_k) + sqrt([sigv^2/dt + (1/12)sigu^2*dt])*Nv
    wtilde_B_BI(:,j+1)    = w_B_BI(:,j+1) + (.5*(bias(:,j+1) + bias(:,j))) + sqrt(((sigv.^2)/(dt)) + (1/12)*(sigu.^2)*dt).*Nv(:,j+1);
 end

Gyro.Measurements  = wtilde_B_BI;
Gyro.Bias          = bias;
Gyro.Nu            = Nu;
Gyro.Nv            = Nv;

end
