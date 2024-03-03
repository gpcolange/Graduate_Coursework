%% Homework 6 Gabriel Colangelo
clear
close all
clc

%% Problem 2
% Polytopic nonlinear system, A(x) = A0 + psi*delta_A
A0          = [-2 1; -1 -3];
delta_A     = [0 1; -1 0];

% Initialize gamma
gamma       = 1;

% Initialize counter
count       = 0;

% Initialize while loop logic
tfeas       = -1;

% Options for feasp - silent
opts        = [0;0;0;0;1];

% Create iterative loop 
while tfeas < 0

    % Increase counter
    count   = count + 1;
   
    % Create counter break
    if count > 1000
        disp('Supremal value of gamma not found')
        fprintf('\n')
        break
    end

    % LMI toolbox setup
    setlmis([]);

    % Matrices for extreme values of psi, a = 0, b = gamma
    A1      = A0;
    A2      = A0 + gamma*delta_A;
    
    % Positive definite matrix
    P       = lmivar(1, [2,1]);

    % Create LMI's: P*Ai + Ai'*P < 0
    lmi1    = newlmi;
    lmiterm([lmi1,1,1,P],1,A1,'s');

    lmi2    = newlmi;
    lmiterm([lmi2,1,1,P],1,A2,'s');

    Plmi    = newlmi;
    lmiterm([-Plmi,1,1,P],1,1);
    lmiterm([Plmi,1,1,0],1);
    lmis    = getlmis;

    % Solve LMIS
    [tfeas, xfeas] = feasp(lmis,opts);

    % Create P matrix
    P       = dec2mat(lmis,xfeas,P);

    % If feasible increase gamma, save latest gamma
    if tfeas < 0
        gamma_max   = gamma;
        gamma       = gamma + .01;
    end

end

% Check P and lyapunov equation
disp('Maximum P is')
disp(P)

disp('Determinant of P is')
disp(det(P))

% Output LMI's
disp('P*A1 + A1''P =')
disp(P*A1 + A1'*P)

disp('P*A2 + A2''P =')
disp(P*A2 + A2'*P)

disp('Eigenvalues of P*A1 + A1''P')
disp(eig(P*A1 + A1'*P))

disp('Eigenvalues of P*A2 + A2''P')
disp(eig(P*A2 + A2'*P))

%% Problem 3
% Polytopic nonlinear system, A(x) = A0 + psi*delta_A
A0          = [0 1; -2 -1];
delta_A     = [0 0; 1 0];

% Set Gamma
gamma       = 1;

% Matrices for extreme values of psi, a = -gamma, b = gamma
A1          = A0 - gamma*delta_A;
A2          = A0 + gamma*delta_A;

% Initialize alpha
alpha       = 0;    

% Initialize counter
count       = 0;

% Initialize while loop logic
tfeas       = -1;

% Create iterative loop 
while tfeas < 0

    % Increase counter
    count   = count + 1;
   
    % Create counter break
    if count > 1000
        disp('Supremal value of alpha not found')
        fprintf('\n')
        break
    end

    % LMI toolbox setup
    setlmis([]);
    
    % Positive definite matrix
    P       = lmivar(1, [2,1]);

    % Create LMI's: P*Ai + Ai'*P <= -2*alpha*P
    lmi1    = newlmi;
    lmiterm([lmi1,1,1,P],1,A1,'s');
    lmiterm([-lmi1 1 1 P], -2*alpha,1); % -2*alpha*P term, RHS

    lmi2    = newlmi;
    lmiterm([lmi2,1,1,P],1,A2,'s');
    lmiterm([-lmi2 1 1 P], -2*alpha,1); % -2*alpha*P term, RHS

    Plmi    = newlmi;
    lmiterm([-Plmi,1,1,P],1,1);
    lmiterm([Plmi,1,1,0],1);

    lmis    = getlmis;

    % Solve LMIS
    [tfeas, xfeas] = feasp(lmis,opts);

    % Create P matrix
    P       = dec2mat(lmis,xfeas,P);

    % If feasible increase alpha, save latest alpha
    if tfeas < 0
        alpha_max   = alpha;
        alpha       = alpha + .001;
    end
end

% Check P and lyapunov equation
disp('Maximum P is')
disp(P)

disp('Eigenvalues of P are')
disp(eig(P))

% Check LMI solver results
disp('Eigenvalues of P*A1 + A1''*P + 2*alpha*P are')
disp(eig(P*A1 + A1'*P + 2*alpha_max*P))

disp('Eigenvalues of P*A2 + A2''*P + 2*alpha*P are')
disp(eig(P*A2 + A2'*P + 2*alpha_max*P))

fprintf('The largest rate of exponential convergence is %.3f \n',alpha_max)

%% Problem 4

% State dependent A(x) = A0 + psi1*delta_A1 + psi2*delta_A2
delta_A1        = zeros(4);
delta_A2        = zeros(4);
delta_A1(3,1)   = 1;
delta_A2(4,2)   = 1;

% Bounds on psi1 and psi2
a1              = -1;
b1              = 1;
a2              = -1;
b2              = 1;

% Initialize spring constant
K               = 0.1;

% Initialize counter
count           = 0;

% Initialize while loop logic
tfeas           = 1;

% Create iterative loop 
while tfeas > 0

    % Increase counter
    count   = count + 1;
   
    % Create counter break
    if count > 1000
        disp('Stable spring constant not found')
        fprintf('\n')
        break
    end

    % LMI toolbox setup
    setlmis([]);

    % Constant matrix
    A0      = [0 0 1 0; 0 0 0 1; -2*K K -2 1; K -K 1 -1];

    % Extreme matrices
    A1      = A0 + a1*delta_A1 + a2*delta_A2;
    A2      = A0 + a1*delta_A1 + b2*delta_A2;
    A3      = A0 + b1*delta_A1 + a2*delta_A2;
    A4      = A0 + b1*delta_A1 + b2*delta_A2;
 
    % Positive definite matrix
    P       = lmivar(1, [4,1]);

    % Create LMI's: P*Ai + Ai'*P < 0
    lmi1    = newlmi;
    lmiterm([lmi1,1,1,P],1,A1,'s');

    lmi2    = newlmi;
    lmiterm([lmi2,1,1,P],1,A2,'s');
     
    lmi3    = newlmi;
    lmiterm([lmi3,1,1,P],1,A3,'s');
     
    lmi4    = newlmi;
    lmiterm([lmi4,1,1,P],1,A4,'s');

    Plmi    = newlmi;
    lmiterm([-Plmi,1,1,P],1,1);
    lmiterm([Plmi,1,1,0],1);
    lmis    = getlmis;

    % Solve LMIS
    [tfeas, xfeas] = feasp(lmis,opts);

    % Create P matrix
    P       = dec2mat(lmis,xfeas,P);

    % If not feasible, increase K
    if tfeas > 0
       K       = K + .05;
    end

end

% Check P and lyapunov equation
fprintf('\n')
disp('Final P is')
disp(P)

disp('Eigenvalues of P are')
disp(eig(P))

% Output LMI's
disp('Eigenvalues of P*A1 + A1''P')
disp(eig(P*A1 + A1'*P))

disp('Eigenvalues of P*A2 + A2''P')
disp(eig(P*A2 + A2'*P))

disp('Eigenvalues of P*A3 + A3''P')
disp(eig(P*A3 + A3'*P))

disp('Eigenvalues of P*A2 + A2''P')
disp(eig(P*A3 + A3'*P))

fprintf(['A spring constant value that guarantees' ...
        ' the system is globally exponentially stable about' ...
        ' the zero solution is K = %.1f \n'],K)
