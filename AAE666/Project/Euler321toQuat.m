function [q] = Euler321toQuat(phi, theta, psi)
% Inputs: (phi, theta, psi)
% Roll  - phi   [rad]
% Pitch - theta [rad]
% Yaw   - psi   [rad]

% Shorthand for cos and sine
c       = @(x) cos(x);
s       = @(x) sin(x);

% Create DCM: Nav to Body
DCM_B_N = [ c(psi)*c(theta), c(theta)*s(psi), -s(theta);...
            c(psi)*s(theta)*s(phi) - s(psi)*c(phi), c(psi)*c(phi) + s(psi)*s(theta)*s(phi), c(theta)*s(phi);...
            s(psi)*s(phi) + c(psi)*s(theta)*c(phi), -c(psi)*s(phi) + s(psi)*s(theta)*c(phi), c(theta)*c(phi)];

% Extract quaternion: desribes attitude matrix for nav to body
q4      = sqrt(1 + trace(DCM_B_N))/2;
q1      = (DCM_B_N(2,3) - DCM_B_N(3,2))/(4*q4);
q2      = (DCM_B_N(3,1) - DCM_B_N(1,3))/(4*q4);
q3      = (DCM_B_N(1,2) - DCM_B_N(2,1))/(4*q4);
q       = [q1;q2;q3;q4];

% Verify quaternion gives attitude matrix (Nav to Body)
A_B_N   = [(q1^2 - q2^2 - q3^2 + q4^2), 2*(q1*q2 + q3*q4), 2*(q1*q3 - q2*q4);...
           2*(q1*q2 - q3*q4), (-q1^2 + q2^2 - q3^2 + q4^2), 2*(q2*q3 + q1*q4);...
           2*(q1*q3 + q2*q4), 2*(q2*q3 - q1*q4), (-q1^2 - q2^2 + q3^2 + q4^2)];

% Check attitude matrix properties
if (max(max(DCM_B_N - A_B_N)) > 1e-6) || (abs(det(A_B_N) - 1) > 1e-6) 
    disp('Attitude Matrix is incorrect')
end

end