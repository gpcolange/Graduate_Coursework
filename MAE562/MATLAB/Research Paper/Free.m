%%% Plant Dynamics
function zdot  = Free(t,z)
z1          = z(1,1); % z1 = r
z2          = z(2,1); % z2 = rdot
z3          = z(3,1); % z3 = thetadot

% Equations of motion is first order form
zdot(1,1)   = z2;
zdot(2,1)   = z1*z3^2;
zdot(3,1)   = z2*z3/z1;
end