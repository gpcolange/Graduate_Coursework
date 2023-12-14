clear
close all
clc

%% Unconstrained MPC
% Linearized Model
A = [0 1 0 0; 5.0449 0 -.9992 0;0 0 0 1; -6.9367 0 0 0];
B = [0 0;.0455 0; 0 0; 0 1/3];
C = [1 0 0 0; 0 0 1 0];
D = [0;0];

m1 = 10;
r1  = 1;
m2 = 3;
g = 9.81;

tf = 10;

% Equilibrium Pair
xe = [pi/4 0 2 0]';
ue = [113.1371;21.2132]';

% step size
h  = .05;

% discrete model
[phi,gamma] = c2d(A,B,h);

% dimensions
[p,n] = size(C);
[~,m] = size(gamma);

% IC
x0 = [deg2rad(-30) 0 1 0]';

% Prediction Horizon
Np = 25;

% target and reference vector
target = [pi/4 2]';
rp     = kron(ones(Np,1),target);

% weights
R      = .05*eye(2*Np);
Q      = 5*eye(2*Np);

% W and Z Calculation
% phi_a   = eye(n+p,n+p);
% gamma_a = zeros(n+p,m);
phi_a  = [phi zeros(4,2); C*phi eye(2,2)];
gamma_a= [gamma;C*gamma];
C_a  = [zeros(2,4), eye(2,2)];

W = [];
for i = 1:Np
    W = [W;C_a*phi_a^i];
end

Z = zeros(Np*p, Np*m);
Z(1:p,1:m) = C_a*gamma_a;
temp = C_a*gamma_a;

for i = 1:Np-1
    temp = [C_a*phi_a^i*gamma_a temp];
    Z(i*p+1:(i+1)*p,1:size(temp,2)) = temp;
end

x = x0;
u = [0 0]';
y = C*x;
out = [y];

t = 0;
time = t;
i = 2;
while t < tf
    t = t + h;

    %% Plant Dynamics
    numf = (-2*m2*x(2)*x(3)*x(4)-g*cos(x(1))*(m1*r1+m2*x(3)) + u(1));
    denf = (m1*r1^2+m2*x(3)^2);
    f2   = numf/denf;
    xdot = [x(2);f2;x(4);x(2)^2*x(3) - g*sin(x(1))+u(2)/m2];
    xnew = x + xdot*h;
    y = C*xnew;
    out = [out y];

    %%
    % Augmented State Vector
    xa = [xnew - x;y];
    x = xnew;
    du = [eye(2) zeros(2,2*Np-2)]*inv(R + Z'*Q*Z)*Z'*Q*(rp - W*xa);
    u  = u + du;
    time = [time t];
    i = i + 1;
end

figure
plot(time,out)

%% ReRun w constraint
x = x0;
u = [0 0]';
y2(:,1) = C*x;

t = 0;
i = 2;

Dg          = [-Z;Z];
dU          = zeros(m*Np,1);
mu          = zeros(length(Dg),1);

Ymax        = kron(ones(Np,1),[10*pi 2]');
Ymin        = kron(ones(Np,1),[-10*pi 1]');


% Optimization Parameters
alpha       = .05;
beta        = alpha;
i = 2;

while t < tf
    t = t + h;

    %% Plant Dynamics
    numf = (-2*m2*x(2)*x(3)*x(4)-g*cos(x(1))*(m1*r1+m2*x(3)) + u(1));
    denf = (m1*r1^2+m2*x(3)^2);
    f2   = numf/denf;
    xdot = [x(2);f2;x(4);x(2)^2*x(3) - g*sin(x(1))+u(2)/m2];
    xnew = x + xdot*h;
    y2(:,i) = C*xnew;

    % Augmented State Vector
    xa = [xnew - x;y2(:,i)];
    
    % Optimizer
    for j = 1:250

        % Cost Function Gradient
        grad_J  = -Z'*Q*(rp - W*xa -Z*dU) + R*dU;

        % Constraint Function
        %g       = Dg*dU - [-Umin + U(:,i); Umax - U(:,i); -Ymin + W*xa; Ymax - W*xa];
        gc      = Dg*dU - [-Ymin + W*xa; Ymax - W*xa];

        % 1st Order Lagrangian Algorithim
        dU      = dU - alpha*(grad_J + Dg'*mu);
        mu      = max(mu + beta*gc,0);
    end

%     % Cost Function to be minimized
%     J           = @(du) 1/2*(rp - W*xa - Z*du)'*Q*(rp - W*xa - Z*du) + 1/2*du'*R*du;
% 
%     % Constraints Ac*X <= Bc
%     Ac          = [-Z;Z];
%     Bc          = [-Ymin + W*xa; Ymax - W*xa];
% 
%     % Optimize
%     dU           = fmincon(J,zeros(m*Np,1),Ac,Bc);

    x  = xnew;
    du = [eye(2) zeros(2,2*Np-2)]*dU;
    u  = u + du;
    i = i + 1;
end

figure
plot(time,out,'--r',time,y2,'-k')
legend('No Constraint','No Constraint','Constraint','Constraint')
