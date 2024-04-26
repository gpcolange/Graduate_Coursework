%% AAE666 Project - Nonlinear Dynamic Inversion Control Law with Euler Angles
clear
close all
clc

%% Setup Plant
syms Ixx Iyy Izz Ixz w1 w2 w3 phi theta psi u1 u2 u3 r1ddot r2ddot r3ddot r1dot...
     r1 r2 r3 v1 v2 v3 Kp Kd r2dot r3dot real;

% Feedforward acceleration
rddot   = [r1ddot;r2ddot;r3ddot];

% Auxiliary input
v       = [v1;v2;v3];

% Moment of Inertia Tensor
I       = [Ixx 0 -Ixz; 0 Iyy 0; -Ixz 0 Izz];

% Angular Velocity Vector
w       = [w1;w2;w3];

% Input Vector
u       = [u1;u2;u3];

% 3-2-1 Euler Angles
eul     = [phi;theta;psi];

% Rotational Kinematic Equation d/dt(eul) = C*w
C       = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);...
           0 cos(phi) -sin(phi);...
           0 sin(phi)*sec(theta) cos(phi)*sec(theta)];

% Euler Rotational Equations of Motion
wdot    = simplify(I\(u - cross(w,I*w)));

% Number of states & inputs
x       = [eul;w];
n       = length(x);
m       = length(u);

%% Non-linear Dynamic Inversion/ Input-Output Feedback Linearization

% Get plant into form of xdot = f(x) + g(x)u, y = h(x)
f       = [C*w;-I\cross(w,I*w)];
g       = [zeros(n - length(I),m);inv(I)];
h       = [phi;theta;psi];

% Lie Derivatives - system has relative order 2
Lgh     = jacobian(h,x)*g
Lfh     = jacobian(h,x)*f
LgLfh   = jacobian(Lfh,x)*g
Lf2h    = jacobian(Lfh,x)*f

% Define F & G
F       = Lf2h;
G       = LgLfh;

% Jacobian of inverse kinematics
gamma   = [(w2*cos(phi)*tan(theta) - w3*sin(phi)*tan(theta)), (w2*sin(phi) + w3*cos(phi))/(cos(theta)^2) 0;...
           (-w2*sin(phi) - w3*cos(phi)), 0, 0;
           (w2*cos(phi) - w3*sin(phi))/cos(theta), (sin(theta)/(cos(theta)^2))*(w2*sin(phi) + w3*cos(phi)), 0];

% Control Law: U = inv(G)*(-F + rddot + v)
U       = simplify(G\(-F + rddot + v));
U_hand  = (I*inv(C))*(-gamma*C*w + C*inv(I)*(cross(w,I*w))...
          + rddot + v);
simplify(U - U_hand);

% Verify error is linear
yddot   = simplify(F + G*U);

% Control Gain Matrix
K       = [Kp Kd 0 0 0 0; 0 0 Kp Kd 0 0; 0 0 0 0 Kp Kd];

% Error vectors
e       = [(r1 - phi); (r1dot - C(1,:)*w);...
           (r2 - theta); (r2dot - C(2,:)*w);...
           (r3 - psi); (r3dot - C(3,:)*w)];

% Closed loop dynamics
xdot    = subs(simplify(f + g*U),[v1;v2;v3],K*e)

