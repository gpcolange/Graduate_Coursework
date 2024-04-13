function [phi_calc,theta_calc, psi_calc] = QuaterniontoEuler321(e1,e2,e3,e4)
% Convert quaternion to 3-2-1 Euler Angles

% Initialize Euler angles
phi_calc    = zeros(length(e1),1);
theta_calc  = phi_calc;
psi_calc    = phi_calc;

% Extract Euler Angles from 3-2-1 Attitude Matrix
for i = 1:length(e1)
    q1              = e1(i);
    q2              = e2(i);
    q3              = e3(i);
    q4              = e4(i);
    A               = [(q1^2 - q2^2 - q3^2 + q4^2), 2*(q1*q2 + q3*q4), 2*(q1*q3 - q2*q4);...
                       2*(q1*q2 - q3*q4), (-q1^2 + q2^2 - q3^2 + q4^2), 2*(q2*q3 + q1*q4);...
                       2*(q1*q3 + q2*q4), 2*(q2*q3 - q1*q4), (-q1^2 - q2^2 + q3^2 + q4^2)];
    psi_calc(i)     = atan(A(1,2)/A(1,1));
    theta_calc(i)   = -asin(A(1,3));
    phi_calc(i)     = atan(A(2,3)/A(3,3));
end

end