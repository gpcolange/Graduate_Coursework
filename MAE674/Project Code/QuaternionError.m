function [delq] = QuaternionError(q,qhat)
% Calculates multiplicative quaternion error between true and estimated
% quaternion. 
% delq      = q x qhat^-1
% where composition of the quaternion is bilinear with
% q' x q    = [E(q) q]q'
% therefore
% delq      = [E(qhat^-1) qhat^-1]q

rhohat      = qhat(1:3);

% qhat^-1 = [-rho,q4]'
qhatinv     = [-rhohat;qhat(4)];

rhohatinv   = qhatinv(1:3);

% E = [q4*I_3x3 + skew(rho); -rho']
E           = [(qhatinv(4)*eye(3)+ skew(rhohatinv));-rhohatinv'];

delq        = [E, qhatinv]*q;

end