%% AAE666 Project - Nonlinear Dynamic Inversion Control Law with quaternions
clear
close all
clc

%% Setup Plant
syms Ixx Iyy Izz Ixz w1 w2 w3 q1 q2 q3 q4 u1 u2 u3 r1ddot r2ddot r3ddot v1 v2 v3 real;

% Moment of Inertia Tensor
I       = [Ixx 0 -Ixz; 0 Iyy 0; -Ixz 0 Izz];

% Angular Velocity Vector
w       = [w1;w2;w3];

% Input Vector
u       = [u1;u2;u3];

% Quaternion vector portion
rho     = [q1;q2;q3];

% Quaternion
q       = [rho;q4];

% Euler Rotational Equations of Motion
wdot    = simplify(I\(u - cross(w,I*w)));

% Define vector dynamics portion of quaternion kinematics
alpha   = q4*eye(3) + skew(rho);

% Quaternion kinematics matrix: qdot = 1/2*E(q)*w
E       = [alpha;-rho'];

% Number of states & inputs
x       = [q;w];
n       = length(x);
m       = length(u);

%% Non-linear Dynamic Inversion/ Input-Output Feedback Linearization

% Get plant into form of xdot = f(x) + g(x)u, y = h(x)
f       = [1/2*E*w;-I\cross(w,I*w)];
g       = [zeros(n - length(I),m);inv(I)];
h       = rho;

% Lie Derivatives - system has relative order 2
Lgh     = jacobian(h,x)*g
Lfh     = jacobian(h,x)*f
LgLfh   = jacobian(Lfh,x)*g
Lf2h    = jacobian(Lfh,x)*f

% Define F & G
F       = Lf2h;
G       = LgLfh;

% Control Law: U = inv(G)*(-F + rddot + v)
rddot   = [r1ddot;r2ddot;r3ddot];
v       = [v1;v2;v3];
U       = simplify(G\(-F + rddot + v));
U_hand  = (2*I*inv(alpha))*(1/4*skew(w)*alpha*w + 1/4*w*rho'*w + ...
          1/2*alpha*inv(I)*(cross(w,I*w)) + rddot + v);
simplify(U - U_hand);

% Verify yddot = rddot + v
yddot   = simplify(F + G*U)

% Zero Dynamics: y = ydot = 0, yddot = 0 = rddot + v
xdot0   = subs(f + g*U,[rho;r1ddot;r2ddot;r3ddot;v1;v2;v3],zeros(9,1))

% Update zero dynamics with q4 ~= 0, w = wdot = 0;
xdot0   = subs(xdot0,w,zeros(3,1))
