clear
close all
clc

%% Problem 1

%%% State Space Matrices
A       = [0 -1; 1 -1];
B       = [1; 0];
C       = eye(2);               % Output both states
D       = zeros(2,1);

t       = (0:.005:2)';          % time vector
x0      = zeros(2,1);           % Initial contiosn x(0) = 0
u       = ones(length(t),1);    % Step Input

sys     = ss(A,B,C,D);          % System state space object

[~,T,X] = lsim(sys,u,t,x0);     % lsim output

%%% Hand Written Solutions for State Trajectories
x1      = exp(-t/2).*(-cos((sqrt(3)*t)/2) + (1/sqrt(3))*sin((sqrt(3)*t)/2)) + 1;
x2      = 1 - exp(-t/2).*(1/sqrt(3) * sin((sqrt(3)*t)/2) + cos((sqrt(3)*t)/2));

figure
subplot(3,1,1)
plot(T,X(:,1),'-r',t,x1,'--k')
xlabel('Time [s]')
title('State 1 Trajectory')
legend('lsim','handwritten')
grid minor
ylabel('x_1')
subplot(3,1,2)
plot(T,X(:,2),'-r',t,x2,'--k')
title('State 2 Trajectory')
xlabel('Time [s]')
grid minor
ylabel('x_2')
subplot(3,1,3)
plot(x1,x2)
title('x2 vs x1')
grid minor
xlabel('x_1')
ylabel('x_2')
sgtitle('Problem 1')


%% Problem 2

Ts      = 0.2;  % step size [s]

disp('Hand Written Discrete State Space Model')

% Hand Written discrete A matrix
Ad      = exp(-.2/2)*[cos(.2*sqrt(3)/2)+ (1/sqrt(3)*sin(.2*sqrt(3)/2)),...
          -2*sin(.2*sqrt(3)/2)/sqrt(3);2*sin(.2*sqrt(3)/2)/sqrt(3),...
          cos(.2*sqrt(3)/2)- (1/sqrt(3)*sin(.2*sqrt(3)/2))]
 
% Hand Written discrete B matrix
Bd      = [exp(-.2/2)*(-cos(.2*sqrt(3)/2) + sin(.2*sqrt(3)/2)/sqrt(3)) + 1;...
          -exp(-.2/2)*(sin(.2*sqrt(3)/2)/sqrt(3) + cos(.2*sqrt(3)/2)) + 1]
      
% Hand Written discrete C matrix
Cd      = [0 1]

% Hand Written discrete D matrix
Dd      = 0

% Continous time ss object
sysc    = ss(A,B,[0 1],0);

disp('c2d Discrete State Space Model')
sysd            = c2d(sysc, Ts, 'zoh')          % c2d discrete ss object

time            = (0:Ts:2)';                    % Time Vector
K               = round(time/Ts);               % Integer Steps
x               = zeros(2,length(K));           % Initialize State Vector
x(:,1)          = x0;
u               = ones(1,length(K));            % Step Input
zeroStateSum    = 0;
y               = zeros(1,length(K));

% x[k] = A^k*x0 + sum(Ad^(k-1-i))*B*u[i] + D*u[k]
for count = 2:length(K)
    k               = K(count);
    zeroStateSum    = Ad^(k - 1)*Bd*u(k) + zeroStateSum;
    x(:,count)      = (Ad^k)*x0 + zeroStateSum;
    y(:,count)      = Cd*x(:,count);
end

% Discrete lsim output
[Y,T,X]     = lsim(sysd,u,time,x0);

figure
subplot(3,1,1)
hold on
stairs(T,X(:,1),'-r')
stairs(time,x(1,:),'--k')
hold off
grid minor
title('State 1 Trajectory, discrete ZOH')
ylabel('x_1')
legend('lsim','hand written')
subplot(3,1,2)
hold on
stairs(T,X(:,2),'-r')
stairs(time,x(2,:),'--k')
hold off
grid minor
title('State 2 Trajectory, discrete ZOH')
ylabel('x_2')
subplot(3,1,3)
hold on
stairs(T,Y,'-r')
stairs(time,y,'--k')
hold off
grid minor
title('Output vs time, discrete ZOH')
ylabel('y')
xlabel('time [s]')
sgtitle('Problem 2')

figure
hold on
step(sysc,2)
step(sysd,2)
grid minor
title('Step Response for Continous and Discrete System')
legend('Continous','Discrete')

%% Problem 5

% Continous Time LTI A Matrix
A_CT        = [-2 0 0; 1 0 1; 0 -2 -2];

% Symmetric Positive Definite Matrix
Q           = eye(size(A_CT));

% Solve Lyapunov Equation for P
P_CT         = lyap(A_CT',Q)

% Check Eigenvalues of P for Positive Definite
lambda_CT    = eig(P_CT)

if (lambda_CT > 0)
    disp('P is positive definite, Continous LTI System is asymptotically stable')
else
    disp('P is not positive definite, Continous LTI System is not asymptotically stable')
end

%% Problem 6

% Discrte Time LTI A matrix
A_DT        = [-.8 0 0; .4 0 .4; 0 -.8 -.8];

% Symmetric Positive Definite Matrix
Q           = eye(size(A_DT));

% Solve Discrete Lyapunov Equation for P
P_DT        = dlyap(A_DT',Q)

% Check Eigenvalues of P for Positive Definite
lambda_DT   = eig(P_DT)

if (lambda_DT > 0)
    disp('P is positive definite, Discrete LTI System is asymptotically stable')
else
    disp('P is not positive definite, Discrete LTI System is not asymptotically stable')
end