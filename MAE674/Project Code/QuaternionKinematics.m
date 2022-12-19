function qdot = QuaternionKinematics(t,q,w)
% Quaternion Kinematic Differential Equation

    % qdot = 1/2*Omega(w)*q
    Omega       = [ -skew(w), w;...
                    -w', 0];
    
    qdot        = (1/2)*Omega*q;

end