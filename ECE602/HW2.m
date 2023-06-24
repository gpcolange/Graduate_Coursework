clear
clc

%% Problem 4 Code

% Matrix of coeffcicients to solve for remainder polynomial 
M       = [4 2 1; 4 -2 1; 4 1 0];  
Minv    =  inv(M)

% Problem 4 A matrix 
A       = [0 1 3; -2 1 1;2 1 1];
[T,J]   = jordan(A)

% Verify A = T*J*T^-1
A_check = T*J*inv(T)

%% Problem 5 Code

% e^A from Laplace Inverse and Jordan Normal Form
eA_hand = [exp(1), 0; 3/2*(exp(1) - exp(-1)), exp(-1)]

% Problem 5 A Matrix
A       = [1 0; 3 -1];

% Matrix exponential
eA_mat  = expm(A)