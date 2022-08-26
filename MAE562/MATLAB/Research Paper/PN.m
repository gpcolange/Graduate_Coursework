%%% Proportional Navigation
function zdot  = PN(t,z,lambda)
z1          = z(1,1); % z1 = r
z2          = z(2,1); % z2 = rdot
z3          = z(3,1); % z3 = theta
z4          = z(4,1); % z4 = thetadot

% Equations of motion is first order form
zdot(1,1)   = z2;
zdot(2,1)   = (1-lambda)*z1*z4^2;
zdot(3,1)   = z4;
zdot(4,1)   = (lambda-2)*z2*z4/z1;
end