clear
close all
warning off
clc

% Symbolic variables
syms x1 x2 x3 x4 x5 x6 u u1 u2 s x(t) x_dot(t) theta1(t) theta_dot_1(t)...
     theta2(t) theta_dot_2(t) t theta_ddot_1 x_ddot theta_ddot_2...
     m1 m2 M L1 L2 g real

%% Problem 1 - Open Loop Lyapunov Stability
disp('-----------------------Begin Problem 1-----------------------------')
fprintf('\n')
% state vector
% x1    = x
% x2    = theta1
% x3    = theta2
% x4    = x_dot
% x5    = theta_dot_1
% x6    = theta_dot_2
x_state = [x1;x2;x3;x4;x5;x6];

% System Parameters
m1_num  = 0.5;
L1_num  = 0.5;
m2_num  = 0.75;
L2_num  = 0.75;
M_num   = 1.5;
g_num   = 9.81;

% Non-linear Model
xdot    = DIPC([],x_state,u,m1_num,m2_num,M_num,L1_num,L2_num,g_num);

% Ouputs
y       = [x1;x2;x3];

% Equilibirum pair - origin
xe       = zeros(6,1);
ue       = 0;

% Jacobian Matrices/ Linearized Model about origin
A       = double(subs(jacobian(xdot,[x1;x2;x3;x4;x5;x6]),[x_state;u],[xe;ue]));
B       = double(subs(jacobian(xdot,u),[x_state;u],[xe;ue]));
C       = double(jacobian(y,[x1;x2;x3;x4;x5;x6]));
D       = double(jacobian(y,u));

% Lyapunov Equation: A'P + PA = -Q
Q       = eye(6);

% Create symbolic symmetric matrix P
P       = sym('P',[6 6],'real');
P       = tril(P,0) + tril(P,-1).';

% Define Symbolic Lyapunov
Lyap_eqn= A'*P + P*A == -Q;

% Get system of equations
eqns    = tril(Lyap_eqn);
eqns    = eqns(eqns~=0);

% Get vector of unknowns
vars    = tril(P);
vars    = vars(vars~=0);

% Solve lyapunov equation
sol     = solve(eqns,vars);

% Extract solutions
fnames  = fieldnames(sol);
for i = 1:length(fnames)
   sol_vec = double(sol.(fnames{i}));
end

% Check if solution to Lyapunov equation is empty & A's eigenvalues are in RHP
if (isempty(sol_vec) == 1) && (max(eig(A) > 0) == 1)
    disp('A solution to the continous time Lyapunov matrix equation does not exist');
    disp('This is because the system is unstable as the eigenvalues of A are in the right hand plane');
    disp('Thus the equilibrium state/open loop system is NOT asymptotically stable in the sense of Lyapunov')
end

%% Problem 2 - Linear State Feedback Controller Design
disp('-----------------------Begin Problem 2-----------------------------')
fprintf('\n')

% % Desired closed loops poles
% pole_cl     = [-1, -1.5, -1.75, -2, -2.5, -3];
% 
% disp('The linear state-feedback controller for the linearized model is: del_u = -K*del_x')
% disp('The control gains for the control law del_u = -K*del_x are:')
% 
% % Control gains
% K           = acker(A,B,pole_cl)
% 
% disp('This gain set will solve the pole placement problem for the linearzied model')

% Get dimensions of B
[n, m] = size(B);

% Use CVX to solve matrix inequality and determine K
cvx_begin sdp quiet

% Variable definition
variable S(n, n) symmetric
variable Z(m, n)

% LMIs
S*A' + A*S -Z'*B'- B*Z +2*S <= -eps*eye(n);
S >= eps*eye(n);
cvx_end

disp('The linear state-feedback controller for the linearized model is: del_u = -K*del_x')
disp('The control gains for the control law del_u = -K*del_x are:')

% compute K matrix
K = Z/S


%% Problem 3 - Closed Loop Transfer Function for Linearized Model
disp('-----------------------Begin Problem 3-----------------------------')
fprintf('\n')

% % Laplace Variable
s           = tf('s');

% Closed Loop transfer Function Matrix Equation
Y_R         = (C - D*K)*inv(s*eye(size(A)) - A + B*K)*B + D;

disp('Transfer function for X to r:')
minreal(Y_R(1),1e-5)

disp('Trasnfer Function for theta_1 to r:')
minreal(Y_R(2),1e-5)

disp('Transfer Function for theta_2 to r:')
minreal(Y_R(3),1e-5)

% disp('Closed loop Transfer function for X to r:')
% simplify(Y_R(1))
% 
% disp('Closed loop Trasnfer Function for theta_1 to r:')
% simplify(Y_R(2))
% 
% disp('Closed loop Transfer Function for theta_2 to r:')
% simplify(Y_R(3))

%% Problem 4 - Closed Loop Lyapunov Function for Linearized Model
disp('-----------------------Begin Problem 4-----------------------------')
fprintf('\n')

% Closed loop A matrix
A_cl    = (A - B*K);

% Lyapunov Function: V = del_x'*P*del_x
disp('The Lyapunov function for the closed-loop system comprised of the linearized model is: V = del_x^T*P*del_x')
disp('Where P is:')

% Solve closed Loop Lyapunov Matrix Equation:A_cl'*P_cl + P_cl*A_cl = -Q
P_cl    = lyap(A_cl',Q);
disp(P_cl)

if min(eig(P_cl) > 0) == 1 && issymmetric(P_cl) == 1
    disp('P is symmetric positive definite')
    disp('Thus the equilibrium state of interest of the closed-loop system is asymptotically stable in the sense of Lyapunov')
end

%% Problem 5 - State Feedback Controller for Two Input System
disp('-----------------------Begin Problem 5-----------------------------')
fprintf('\n')

% Call Lagrangian for DIPC 
L = DIPC_Lagrangian(t,x,x_dot, theta1, theta_dot_1, theta2, theta_dot_2, M, m1,m2, L1, L2, g);

% Solve Lagrange's Equations of Motion
% q = x, Q = u1
eqn_x       = subs(simplify(diff(diff(L,x_dot),t) - diff(L,x)),[diff(x(t),t), diff(theta1(t),t)...
             ,diff(theta2(t),t), diff(x_dot,t), diff(theta_dot_1(t), t), diff(theta_dot_2(t), t)],...
             [x_dot, theta_dot_1, theta_dot_2, x_ddot, theta_ddot_1, theta_ddot_2]) == u1;

% q = theta1, Q = u2
eqn_theta1  = subs(simplify(diff(diff(L,theta_dot_1),t) - diff(L,theta1)),[diff(x(t),t), diff(theta1(t),t)...
             ,diff(theta2(t),t), diff(x_dot,t), diff(theta_dot_1(t), t), diff(theta_dot_2(t), t)],...
             [x_dot, theta_dot_1, theta_dot_2, x_ddot, theta_ddot_1, theta_ddot_2]) == u2;

% q = theta2, Q = 0
eqn_theta2  = subs(simplify(diff(diff(L,theta_dot_2),t) - diff(L,theta2)),[diff(x(t),t), diff(theta1(t),t)...
             ,diff(theta2(t),t), diff(x_dot,t), diff(theta_dot_1(t), t), diff(theta_dot_2(t), t)],...
             [x_dot, theta_dot_1, theta_dot_2, x_ddot, theta_ddot_1, theta_ddot_2]) == 0;

% Solve system of equations for 2nd derivative of states
sys_eqn     = solve([eqn_x,eqn_theta1,eqn_theta2],[x_ddot,theta_ddot_1,theta_ddot_2]);

% Put EOM into state space form
x1_dot      = x4;
x2_dot      = x5;
x3_dot      = x6;
x4_dot      = subs(simplify(sys_eqn.x_ddot),[x theta1 theta2 x_dot theta_dot_1 theta_dot_2],[x1 x2 x3 x4 x5 x6]);
x5_dot      = subs(simplify(sys_eqn.theta_ddot_1),[x theta1 theta2 x_dot theta_dot_1 theta_dot_2],[x1 x2 x3 x4 x5 x6]);
x6_dot      = subs(simplify(sys_eqn.theta_ddot_2),[x theta1 theta2 x_dot theta_dot_1 theta_dot_2],[x1 x2 x3 x4 x5 x6]);

% Define Non-linear system
f           = [x1_dot;x2_dot;x3_dot;x4_dot;x5_dot;x6_dot];
h           = [x1;x2;x3];

% Jacobian Matrices
df_dx       = subs(jacobian(f,[x1;x2;x3;x4;x5;x6]),[m1 m2 M L1 L2 g],[m1_num, m2_num,M_num,L1_num, L2_num, g_num]);
df_du       = subs(jacobian(f,[u1;u2]),[m1 m2 M L1 L2 g],[m1_num, m2_num,M_num,L1_num, L2_num, g_num]);

% Input for equilibirum at origin
u1e         = 0;
u2e         = 0;

disp('The updated linearized model with two inputs is:')

% Redefine State Matrices with extra input
A           = double(subs(df_dx,[x1;x2;x3;x4;x5;x6;u1;u2],[xe;u1e;u2e]))
B           = double(subs(df_du,[x1;x2;x3;x4;x5;x6;u1;u2],[xe;u1e;u2e]))
C
D           = double(jacobian(h,[u1;u2]))

% Get dimensions of new B
[n, m] = size(B);

% Use CVX to solve matrix inequality and determine new K
cvx_begin sdp quiet

% Variable definition
variable S(n, n) symmetric
variable Z(m, n)

% LMIs
S*A' + A*S -Z'*B'- B*Z +2*S <= -eps*eye(n);
S >= eps*eye(n);
cvx_end

disp('The linear state-feedback controller for the new linearized model is: del_u = -K*del_x')
disp('The new control gains for the control law del_u = -K*del_x are:')

% compute K matrix
K = Z/S

%% Problem 6 - State Observer Design
disp('-----------------------Begin Problem 6-----------------------------')
fprintf('\n')

[p, n] = size(C);

% Use CVX to solve matrix inequality and determine L 
cvx_begin sdp quiet

% Variable definition
variable P(n, n) symmetric
variable Y(n, p)

% LMI with robustness term (all eigenvalues less than -2)
A'*P + P*A - C'*Y' - Y*C + 4*P <= -eps*eye(n);
P >= eps*eye(n) 
cvx_end

% solver for observer gain matrix
L  = P\Y

%% Problem 7 - Lyapunov Function for combined observer controller compensator
disp('-----------------------Begin Problem 7-----------------------------')
fprintf('\n')

% A matrix for closed loop system driven by the combined observer controller compensator
A_cl_full = [A, -B*K; L*C, A - L*C - B*K];

% Solve Lyapunov Matrix Equation for combined observer controller compensator system :A_cl'*P_cl + P_cl*A_cl = -Q
P_cl_full = lyap(A_cl_full',eye(12));
disp(P_cl_full)

if min(eig(P_cl_full) > 0) == 1 && issymmetric(P_cl_full) == 1
    disp('P is symmetric positive definite')
    disp('Thus the equilibrium state of interest of the closed-loop system is asymptotically stable in the sense of Lyapunov')
end

%% Problem 8 - Transfer Function for combined observer controller compensator
disp('-----------------------Begin Problem 8-----------------------------')
fprintf('\n')

% Closed Loop transfer Function Matrix Equation
Y_R         = (C - D*K)*inv(s*eye(size(A)) - A + B*K)*B + D;

disp('Obersver - Controller Closed loop Transfer function for X to r:')
minreal(Y_R(1,1))

disp('Obersver - Controller Closed loop Trasnfer Function for theta_1 to r:')
minreal(Y_R(2,1))

disp('Obersver - Controller Closed loop Transfer Function for theta_2 to r:')
minreal(Y_R(3,1))

disp('Obersver - Controller Closed loop Transfer function for X to r2:')
minreal(Y_R(1,2))

disp('Obersver - Controller Closed loop Trasnfer Function for theta_1 to r2:')
minreal(Y_R(2,2))

disp('Obersver - Controller Closed loop Transfer Function for theta_2 to r2:')
minreal(Y_R(3,2))

%% Problem 9 - Simulation

% Observer IC 
z0          = zeros(6,1);

% State IC  [m, rad, rad, m/s, rad/s, rad/s]
x0          = [0 .01 .02 0 0 0]'; 

% Time interval and vector
dt          = 1/200;
time        = (0:dt:5)';

% equilibirum ouput for observer
ye          = subs(h,[x1; x2; x3],xe(1:3));

% ODE45 solver options
options     = odeset('AbsTol',1e-8,'RelTol',1e-8);

% ODE45 Function call
[~, X] = ode45(@(t,x) ControlledDIPC(t, x, A, B, C, D, K, L, [u1e;u2e], ye, M_num, m1_num ,m2_num, L1_num, L2_num, g_num),...
         time,[x0;z0], options);

figure
plot(time, X(:,1:3))
grid on

% figure
% subplot(611)
% plot(time,X(:,1),time,X(:,7))
% subplot(612)
% plot(time,X(:,2),time,X(:,8))
% subplot(613)
% plot(time,X(:,3),time,X(:,9))
% subplot(614)
% plot(time,X(:,4),time,X(:,10))
% subplot(615)
% plot(time,X(:,5),time,X(:,11))
% subplot(616)
% plot(time,X(:,6),time,X(:,12))


%% Problem 9 Part - 2 Animation

% Cart width and height
w           = 1;
h           = .5;

% Graphics handle - cart
cart        = rectangle('position',[X(1,1) - w/2, -h, w, h]);

% Graphics handle - hinge 
hinge      = line('xdata', X(1,1),'ydata',0,'marker','o','markersize',7);

% Graphics handle - mass 1
mass1       = line('xdata', X(1,1) + L1_num*sin(X(1,2)), 'ydata', L1_num*cos(X(1,2)),...
              'marker','o','markersize',10,'MarkerFaceColor','k');

% Graphics handle - bar 1
bar1        = line('xdata', [X(1,1) X(1,1) + L1_num*sin(X(1,2))],'ydata',...
               [0 L1_num*cos(X(1,2))],'linewidth',3);   

% Graphics handle - mass 2
mass2       = line('xdata', X(1,1) + L1_num*sin(X(1,2)) + L2_num*sin(X(1,3)), 'ydata',...
               L1_num*cos(X(1,2))+L2_num*cos(X(1,3)),'marker','o','markersize',10,'MarkerFaceColor','k');

% Graphics handle - bar 2
bar2        = line('xdata', [(X(1,1) + L1_num*sin(X(1,2))), (X(1,1) + L1_num*sin(X(1,2)) + L2_num*sin(X(1,3)))],'ydata',...
               [(L1_num*cos(X(1,2))) (L1_num*cos(X(1,2))+L2_num*cos(X(1,3)))],'linewidth',3);  

h_txt       = text(-1.1,1.1,strcat(['Time = ',' ', num2str(time(1)), ' [s]']));

% Define axis limits
axis([-1.1*(L1_num + L2_num), 1.1*(L1_num + L2_num),-1.1*(L1_num + L2_num), 1.1*(L1_num + L2_num)]);
grid on
xlabel('X [m]')
ylabel('Y [m]')
title('Controlled Double Inverted Pendulum')

% Video stuff
vidobj      = VideoWriter('DIPC.avi');
open(vidobj);
nframes     = length(X);
frames      = moviein(nframes);

for i = 2:nframes

    % Update handles
    set(cart,'position',[X(i,1) - w/2, -h, w, h]);
    set(hinge,'xdata', X(i,1),'ydata',0,'marker','o','markersize',7);
    set(mass1,'xdata', X(i,1) + L1_num*sin(X(i,2)), 'ydata', L1_num*cos(X(i,2)),...
                  'marker','o','markersize',10,'MarkerFaceColor','k');
    set(bar1,'xdata', [X(i,1) X(i,1) + L1_num*sin(X(i,2))],'ydata',...
                   [0 L1_num*cos(X(i,2))],'linewidth',3);   
    set(mass2,'xdata', X(i,1) + L1_num*sin(X(i,2)) + L2_num*sin(X(i,3)), 'ydata',...
                   L1_num*cos(X(i,2))+L2_num*cos(X(i,3)),'marker','o','markersize',10,'MarkerFaceColor','k');
    set(bar2,'xdata', [(X(i,1) + L1_num*sin(X(i,2))), (X(i,1) + L1_num*sin(X(i,2)) + L2_num*sin(X(i,3)))],'ydata',...
                   [(L1_num*cos(X(i,2))) (L1_num*cos(X(i,2))+L2_num*cos(X(i,3)))],'linewidth',3); 
    set(h_txt,'String',strcat(['Time = ',' ', num2str(time(i)), ' [s]']));

    drawnow;
    frames(:,i) = getframe(gcf);
    writeVideo(vidobj,frames(:,i));
end

close(vidobj);

function xdot = DIPC(t,x,u,m1,m2,M,L1,L2,g)

% States and inputs
x1      = x(1,1); % x
x2      = x(2,1); % theta_1
x3      = x(3,1); % theta_2
x4      = x(4,1); % xdot
x5      = x(5,1); % theta_1_dot
x6      = x(6,1); % theta_2_dot


% Equations of Motion
x1dot   = x4;   % xdot
x2dot   = x5;   % theta_1_dot
x3dot   = x6;   % theta_2_dot

% x_ddot
x4dot   = (2*m1*u + m2*u - m2*u*cos(2*x2 - 2*x3) - g*m1^2*sin(2*x2) +...
          2*L1*m1^2*x5^2*sin(x2) - g*m1*m2*sin(2*x2) +...
          2*L1*m1*m2*x5^2*sin(x2) + L2*m1*m2*x6^2*sin(x3) +...
          L2*m1*m2*x6^2*sin(2*x2 - x3))/(2*M*m1 + M*m2 + m1*m2 -...
          m1^2*cos(2*x2) + m1^2 - m1*m2*cos(2*x2) - M*m2*cos(2*x2 - 2*x3));

% theta_1_ddot
x5dot   = -(m1*u*cos(x2) + (m2*u*cos(x2))/2 - (m2*u*cos(x2 - 2*x3))/2 -...
           g*m1^2*sin(x2) - M*g*m1*sin(x2) - (M*g*m2*sin(x2))/2 -...
           g*m1*m2*sin(x2) + (L1*m1^2*x5^2*sin(2*x2))/2 - (M*g*m2*...
           sin(x2 - 2*x3))/2 + (L2*m1*m2*x6^2*sin(x2 + x3))/2 + ...
           L2*M*m2*x6^2*sin(x2 - x3) + (L2*m1*m2*x6^2*sin(x2 - x3))/2 +...
           (L1*m1*m2*x5^2*sin(2*x2))/2 + (L1*M*m2*x5^2*sin(2*x2 - 2*x3))/2)...
           /(L1*(M*m1 + (M*m2)/2 + (m1*m2)/2 - (m1^2*cos(2*x2))/2 + m1^2/2 ...
           - (m1*m2*cos(2*x2))/2 - (M*m2*cos(2*x2 - 2*x3))/2));

% theta_2_ddot
x6dot   = ((m1*u*cos(2*x2 - x3))/2 - (m2*u*cos(x3))/2 - (m1*u*cos(x3))/2 +...
          (m2*u*cos(2*x2 - x3))/2 - (M*g*m1*sin(2*x2 - x3))/2 - ...
          (M*g*m2*sin(2*x2 - x3))/2 + (M*g*m1*sin(x3))/2 + ...
          (M*g*m2*sin(x3))/2 + L1*M*m1*x5^2*sin(x2 - x3) +...
          L1*M*m2*x5^2*sin(x2 - x3) + (L2*M*m2*x6^2*sin(2*x2 - 2*x3))/2)/...
          (L2*(M*m1 + (M*m2)/2 + (m1*m2)/2 - (m1^2*cos(2*x2))/2 + m1^2/2 -...
          (m1*m2*cos(2*x2))/2 - (M*m2*cos(2*x2 - 2*x3))/2));

xdot    = [x1dot;x2dot;x3dot;x4dot;x5dot;x6dot];

end

function L = DIPC_Lagrangian(t,x, x_dot, theta1, theta_dot_1, theta2, theta_dot_2, M, m1,m2, L1, L2, g)

% Lagrangian for DIPC from HW1
L    = (m2*(x_dot(t) + L1*cos(theta1(t))*theta_dot_1(t) + L2*cos(theta2(t))*...
       theta_dot_2(t))^2)/2 + (m1*(x_dot(t) + L1*cos(theta1(t))*...
       theta_dot_1(t))^2)/2 + (m2*(L1*sin(theta1(t))*theta_dot_1(t) +...
       L2*sin(theta2(t))*theta_dot_2(t))^2)/2 + (M*x_dot(t)^2)/2 + ...
       (L1^2*m1*sin(theta1(t))^2*theta_dot_1(t)^2)/2 - ...
       L1*g*m1*cos(theta1(t)) - L1*g*m2*cos(theta1(t)) - L2*g*m2*cos(theta2(t));

end

function xdot = ControlledDIPC(t, x, A, B, C, D, K, L, ue, ye, M, m1,m2, L1, L2, g)
% Define State Vector
x1          = x(1,1);  % x
x2          = x(2,1);  % theta_1
x3          = x(3,1);  % theta_2
x4          = x(4,1);  % xdot
x5          = x(5,1);  % theta_1_dot
x6          = x(6,1);  % theta_2_dot
z1          = x(7,1);  % delta_x1_tilde - estimate of change in x
z2          = x(8,1);  % delta_x2_tilde - estimate of change in theta_1
z3          = x(9,1);  % delta_x3_tilde - estimate of change in theta_2
z4          = x(10,1); % delta_x4_tilde - estimate of change in xdot
z5          = x(11,1); % delta_x5_tilde - estimate of change in theta_1_dot
z6          = x(12,1); % delta_x6_tilde - estimate of change in tbeta_2_dot

% delta_tilde_x - vector of state pertubation estimates
z           = [z1;z2;z3;z4;z5;z6];

% Output vector - x, theta_1, theta_2
y           = [x1;x2;x3];

% Output pertubation vector
del_y       = y - ye;

% Control law
del_u       = - K*z;
u           = del_u + ue;
u1          = u(1);
u2          = u(2);

% State Dynamics
x1dot       = x4;   % xdot
x2dot       = x5;   % theta_1_dot
x3dot       = x6;   % theta_2_dot

x4dot       = ((m2*u2*cos(x2 - 2*x3))/2 - (m2*u2*cos(x2))/2 - m1*u2*cos(x2) +...
           L1*m1*u1 + (L1*m2*u1)/2 - (L1*g*m1^2*sin(2*x2))/2 + L1^2*m1^2*x5^2*sin(x2)...
           - (L1*m2*u1*cos(2*x2 - 2*x3))/2 - (L1*g*m1*m2*sin(2*x2))/2 +...
           L1^2*m1*m2*x5^2*sin(x2) + (L1*L2*m1*m2*x6^2*sin(2*x2 - x3))/2 + ...
           (L1*L2*m1*m2*x6^2*sin(x3))/2)/(L1*(M*m1 + (M*m2)/2 + (m1*m2)/2 -...
           (m1^2*cos(2*x2))/2 + m1^2/2 - (m1*m2*cos(2*x2))/2 - (M*m2*cos(2*x2 - 2*x3))/2));

x5dot       = (M*u2 + m1*u2 + (m2*u2)/2 - (m2*u2*cos(2*x3))/2 - L1*m1*u1*cos(x2) -...
           (L1*m2*u1*cos(x2))/2 + (L1*m2*u1*cos(x2 - 2*x3))/2 + L1*g*m1^2*sin(x2) -...
           (L1^2*m1^2*x5^2*sin(2*x2))/2 - (L1^2*m1*m2*x5^2*sin(2*x2))/2 -...
           (L1^2*M*m2*x5^2*sin(2*x2 - 2*x3))/2 + L1*M*g*m1*sin(x2) +...
           (L1*M*g*m2*sin(x2))/2 + L1*g*m1*m2*sin(x2) +...
           (L1*M*g*m2*sin(x2 - 2*x3))/2 - (L1*L2*m1*m2*x6^2*sin(x2 + x3))/2 -...
           L1*L2*M*m2*x6^2*sin(x2 - x3) - (L1*L2*m1*m2*x6^2*sin(x2 - x3))/2)...
           /(L1^2*(M*m1 + (M*m2)/2 + (m1*m2)/2 - (m1^2*cos(2*x2))/2 + m1^2/2 -...
           (m1*m2*cos(2*x2))/2 - (M*m2*cos(2*x2 - 2*x3))/2));

x6dot       = (m1*u2*cos(x2 + x3) - m1*u2*cos(x2 - x3) - m2*u2*cos(x2 - x3) -...
           2*M*u2*cos(x2 - x3) + m2*u2*cos(x2 + x3) - L1*m1*u1*cos(x3) - ...
           L1*m2*u1*cos(x3) + L1*m1*u1*cos(2*x2 - x3) + L1*m2*u1*cos(2*x2 - x3) -...
           L1*M*g*m1*sin(2*x2 - x3) - L1*M*g*m2*sin(2*x2 - x3) + L1*M*g*m1*sin(x3) +...
           L1*M*g*m2*sin(x3) + 2*L1^2*M*m1*x5^2*sin(x2 - x3) +...
           2*L1^2*M*m2*x5^2*sin(x2 - x3) + L1*L2*M*m2*x6^2*sin(2*x2 - 2*x3))/...
           (L1*L2*(2*M*m1 + M*m2 + m1*m2 - m1^2*cos(2*x2) + m1^2 ...
           - m1*m2*cos(2*x2) - M*m2*cos(2*x2 - 2*x3)));

xdot(1:6,1) = [x1dot;x2dot;x3dot;x4dot;x5dot;x6dot];

% Observer Dynamics
del_y_tilde = C*z + D*del_u;
zdot        = A*z + B*del_u + L*(del_y - del_y_tilde);

xdot(7:12,1)= zdot;
end