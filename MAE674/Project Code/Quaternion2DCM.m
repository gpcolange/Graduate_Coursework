function A = Quaternion2DCM(Quaternion)
%% Converts Quaternion to Attitude Matrix
%% Input must be a 4x1 vector in form of [q1 q2 q3 q4]' where q4 is scalar

% p = [q1 q2 q3]'
rho             = Quaternion(1:3);

% E = [q4*I + [px];-p']
E               = [(Quaternion(4)*eye(3) + skew(rho)); - rho' ]; 

% Psi = [q4*I - [px];-p']
Psi             = [(Quaternion(4)*eye(3) - skew(rho)); - rho' ]; 

% Attitude Matrix
A               = E'*Psi;

end