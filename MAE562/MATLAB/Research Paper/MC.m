%%% Motion Camouflage Navigation
function zdot  = MC(t,z,mu)
% rP          = z(1,1);
% xP          = z(2,1); 
% yP          = z(3,1);
% rT          = z(4,1);
% xT          = z(5,1);
% yT          = z(6,1);
% 
% r           = rp - rT;
% N           = mu*r0;
% lambda      = -(1/norm(r))*(r/norm(r)
% 
% zdot(1,1)   = xP;           % rPdot
% zdot(2,1)   = yP*uP;        % xPdot
% zdot(3,1)   = -xP*up;       % yPdot
% zdot(4,1)   = v*xT;         % rTdot
% zdot(5,1)   = v*yT*uT;      % xTdot
% zdot(6,1)   = -v*xT*uT;     % yTdot

x          = z(1,1);
y          = z(2,1);
xdot       = z(3,1);
ydot       = z(4,1);

% xT          = z(5,1);
% yT          = z(6,1);
% xTdot       = z(7,1);
% yTdot       = z(8,1);

r           = [x; y];
rdot        = [xdot; ydot];

rddot       = mu*((dot((r/norm(r)),rdot))*rdot - (norm(r)^2)*(r/norm(r)));

zdot(1,1)   = xdot;
zdot(2,1)   = ydot;
zdot(3,1)   = rddot(1,1);
zdot(4,1)   = rddot(2,1);

end