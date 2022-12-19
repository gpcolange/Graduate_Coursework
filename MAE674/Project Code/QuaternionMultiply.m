function [delq] = QuaternionMultiply(q1,q2)
% Performs Quaternion Multiplication
% where composition of the quaternion is bilinear with
% q1 x q2    = [E(q2) q2]q1

rhohat      = q2(1:3);

% E = [q4*I_3x3 + skew(rho); -rho']
E           = [(q2(4)*eye(3)+ skew(rhohat));-rhohat'];

delq        = [E, q2]*q1;

end