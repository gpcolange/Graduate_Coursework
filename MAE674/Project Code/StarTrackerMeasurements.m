function StarTracker = StarTrackerMeasurements(TrueQuaternion, R, ReferenceVectors)
%%% Inputs are the true quaternion, Measurement error covariance matrox, and matrix
%%% containing the true refernece vectors of size 3Nxm where N is the number
%%% of vector sensors

Det     = zeros(1,length(ReferenceVectors));
btilde  = zeros(size(ReferenceVectors));
btrue   = btilde;

for i = 1:length(ReferenceVectors)
    A               =  Quaternion2DCM(TrueQuaternion(:,i));
    Det(i)          = det(A);
    
    % ith refence vector
    r               = ReferenceVectors(:,i);
    
    % Noise
    v               = sqrt(diag(R)).*randn(3,1);
    
    % Measurement
    btilde(:,i)     = A*r + v;

    btrue(:,i)      = A*r;
end

StarTracker.ReferenceVectors       = ReferenceVectors;
StarTracker.BodyMeasurements       = btilde;
StarTracker.Determinate            = Det;
StarTracker.BodyTrue               = btrue;

end