function [qc,qcdot, qcddot, wc] = QuaternionCommand(phic,thetac,psic)
%%% Inputs are vectors containing euler angle attitude, velocity, and accelerations
%%% Ex: phic(1) = phi, phic(2) = phidot, phic(3) = phiddot

% Convert Euler angle to quaternion
qc      = Euler321toQuat(phic(1),thetac(1),psic(1));

% Extract euler angles & rates, used in wc & wcdot calculation
theta   = thetac(1);
phi     = phic(1);
phidot  = phic(2);
thetadot= thetac(2);

% Convert Euler rates to angular velocity: wc = C*adot
wc      = [-sin(theta), 0, 1;...
            cos(theta)*sin(phi), cos(phi), 0;...
            cos(phi)*cos(theta), -sin(phi), 0]*[psic(2);thetac(2);phic(2)];

% Quaternion command rate
qcdot   = QuaternionKinematics([],qc,wc);

% Angular acceleration from euler rates and acceleration: 
% wcdot = C*addot + Cdot*adot
wdotc   =   [-sin(theta), 0, 1;...
            cos(theta)*sin(phi), cos(phi), 0;...
            cos(phi)*cos(theta), -sin(phi), 0]*[psic(3);thetac(3);phic(3)]...
           +... 
            [-thetadot*cos(theta), 0 ,0;...
            phidot*cos(phi)*cos(theta) - thetadot*sin(phi)*sin(theta), -phidot*sin(phi), 0;...
            -phidot*sin(phi)*cos(theta) - thetadot*cos(phi)*sin(theta), -phidot*cos(phi),0]...
            *[psic(2);thetac(2);phic(2)];

% Form  rho & E(q) 
% rho = [q1 q2 q3]'
rho     = qc(1:3,1);

% E = [q4*I + [px];-p']
E       = [(qc(4,1)*eye(3) + skew(rho)); - rho' ];

% Quaternion command acceleration
qcddot  = 1/2*E*wdotc + QuaternionKinematics([],qcdot,wc);

end