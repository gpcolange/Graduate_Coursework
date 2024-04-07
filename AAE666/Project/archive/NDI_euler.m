%% AAE666 Project - Nonlinear Dynamic Inversion Control Law with Euler Angles
clear
close all
clc

%% Setup Plant
syms Ixx Iyy Izz Ixz w1 w2 w3 phi theta psi u1 u2 u3 rddot v real;

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

% Control Law: U = inv(G)*(-F + rddot + v)
U       = simplify(G\(-F + rddot + v));
U_hand  = inv(C*inv(I))*(-jacobian(C*w,eul)*C*w + C*inv(I)*(cross(w,I*w))...
          + rddot + v);
simplify(U - U_hand);

