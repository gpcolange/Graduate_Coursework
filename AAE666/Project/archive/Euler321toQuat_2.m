function [q1, q2, q3, q4] = Euler321toQuat_2(phi, theta, psi)
% Inputs: (phi, theta, psi)
% Roll  - phi   [rad]
% Pitch - theta [rad]
% Yaw   - psi   [rad]
% Contains checks for quaternion and attitude definition

% Shorthand for cos and sine
c       = @(x) cos(x);
s       = @(x) sin(x);

% Create DCM: Nav to Body
DCM_B_N = [ c(psi)*c(theta), c(theta)*s(psi), -s(theta);...
            c(psi)*s(theta)*s(phi) - s(psi)*c(phi), c(psi)*c(phi) + s(psi)*s(theta)*s(phi), c(theta)*s(phi);...
            s(psi)*s(phi) + c(psi)*s(theta)*c(phi), -c(psi)*s(phi) + s(psi)*s(theta)*c(phi), c(theta)*c(phi)];

% Create DCM: Body to Nav
C3      = [c(psi), -s(psi), 0; s(psi), c(psi), 0; 0 0 1];
C2      = [c(theta), 0 s(theta); 0 1 0; -s(theta) 0 c(theta)];
C1      = [1 0 0; 0 c(phi) -s(phi); 0 s(phi) c(phi)];
DCM_N_B = C3*C2*C1;

% Check properties of DCM det = 1 and A_B_N' =  A_N_B
if (abs(det(DCM_N_B) - 1) > 1e-6) || (max(max(DCM_N_B ~= DCM_B_N')) == 1)
    disp('DCM properties not satisifed')
end

% Extract quaternion: desribes attitude matrix for nav to body
q4      = sqrt(1 + trace(DCM_B_N))/2;
q1      = (DCM_B_N(2,3) - DCM_B_N(3,2))/(4*q4);
q2      = (DCM_B_N(3,1) - DCM_B_N(1,3))/(4*q4);
q3      = (DCM_B_N(1,2) - DCM_B_N(2,1))/(4*q4);

% Verify quaternion gives attitude matrix (Nav to Body)
A_B_N   = [(q1^2 - q2^2 - q3^2 + q4^2), 2*(q1*q2 + q3*q4), 2*(q1*q3 - q2*q4);...
           2*(q1*q2 - q3*q4), (-q1^2 + q2^2 - q3^2 + q4^2), 2*(q2*q3 + q1*q4);...
           2*(q1*q3 + q2*q4), 2*(q2*q3 - q1*q4), (-q1^2 - q2^2 + q3^2 + q4^2)];
if (max(max(DCM_B_N - A_B_N)) > 1e-6)
    disp('Attitude Matrix is incorrect')
end