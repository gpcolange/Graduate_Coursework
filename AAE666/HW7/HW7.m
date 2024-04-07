%% Homework 7 Gabriel Colangelo
clear
close all
clc

%% Problem 2

% Scaling parameters
lambda1 = 2;
lambda2 = 3;
mu1     = lambda1^2;
mu2     = lambda2^2;

% xdot  = Ax + B1*phi1 + B2*phi2
% zi    = Ci*X
C1      = [1 0 0 0];
C2      = [0 1 0 0];
B1      = [0 0 1 0]';
B2      = [0 0 0 1]';

% Initialize spring constant
K       = 0.1;

% Initialize counter
count   = 0;

% Initialize while loop logic
tfeas   = 1;

% Options for feasp - silent
opts    = [0;0;0;0;1];

% Create iterative loop 
while tfeas > -.01

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
    A       = [0 0 1 0; 0 0 0 1; -2*K K -2 1; K -K 1 -1];
 
    % Positive definite matrix
    P       = lmivar(1, [4,1]);

    % Create LMI -[PA+A'P+C1'C1+C2'C2 PB1 PB2; B1'P -I 0; B2'P 0 -I];
    lmi1    = newlmi;
    lmiterm([lmi1,1,1,P],1,A,'s');      % PA + A'P
    lmiterm([lmi1 1 1 0],mu1*(C1'*C1)); %mu1*C1'C1
    lmiterm([lmi1 1 1 0],mu2*(C2'*C2)); %mu2*C2'C2
    lmiterm([lmi1 1 2 P],1,B1);         %PB1
    lmiterm([lmi1 1 3 P],1,B2);         %PB2
    lmiterm([lmi1 2 2 0],-mu1);         %-mu1*I
    lmiterm([lmi1 3 3 0],-mu2);         %-mu2*I

    Plmi    = newlmi;
    lmiterm([-Plmi,1,1,P],1,1);
    lmiterm([Plmi,1,1,0],1);
    lmis    = getlmis;

    % Solve LMIS
    [tfeas, xfeas] = feasp(lmis,opts);

    % Create P matrix
    P       = dec2mat(lmis,xfeas,P);

    % If not feasible, increase K
    if tfeas > -.01
       K       = K + .1;
    end

end

% Check P
fprintf('\n')
disp('Final P is')
disp(P)

disp('Eigenvalues of P are')
disp(eig(P))

fprintf(['A spring constant value that guarantees' ...
        ' the system is globally exponentially stable about' ...
        ' the zero solution is K = %.1f \n'],K)

% Check QMI
Ctilde = [lambda1*C1;lambda2*C2];
Btilde = [lambda1^-1*B1, lambda2^-1*B2];
Q      = P*A + A'*P + P*Btilde*Btilde'*P + Ctilde'*Ctilde;

fprintf('\n')
disp('The QMI from Equation 14.37 is:')
disp(Q)
disp('With eigenvalues:')
disp(eig(Q))


%% Problem 4

% Laplace variable
s       = tf('s');

% Initial beta 
beta    = 0;

% Initialize counter
count   = 0;

% Initialize Loop logic
logic   = 0;

% Define small posiitve alpha
alpha   = 1e-3;

while logic == 0

    % Increase counter
    count   = count + 1;
   
    % Create counter break
    if count > 1000
        disp('SPR beta not found')
        fprintf('\n')
        break
    end

    % Transfer Funtion
    g       = (beta*s + 1)/(s^2 + s + 2);

    % Transfer function to state space
    [A,B,C,D] = tf2ss(g.Numerator{1},g.Denominator{1});

    % Get size of A
    n        = length(A);

    cvx_begin sdp quiet
    
    % Variable definition
    variable P(n, n) symmetric
    
    % LMIs
    [(P*A + A'*P + 2*alpha*P), (P*B - C');...
     (B'*P - C), -(D + D')] <= -eps*eye(n+1);
    P >= eps*eye(n);
    cvx_end

    % Check in P = P' > 0
    logic = all(eig(P) > 0);

    % If logic is true, increase beta
    if logic == 0
       beta  = beta + .01;
    end
end

% Create controllability and observability gramian
Wc      = gram(ss(A,B,C,D),'c');
Wo      = gram(ss(A,B,C,D),'o');

% Check if gramian is invertible and A is Hurwitz
if det(Wc) ~= 0 && all(eig(A) < 0)
    disp('System is controllable')
end

if det(Wo)~= 0 && all(eig(A) < 0)
    disp('System is observable')
end

fprintf('The minimum beta needed for P = P'' > 0 is %.2f \n',beta)

LMI     = [(P*A + A'*P + 2*alpha*P), (P*B - C');...
           (B'*P - C), -(D + D')]

disp('The eigenvalues of the LMI are: ')
disp(eig(LMI))
