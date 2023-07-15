clear
clc

%% Problem 1
disp('--Problem 1 --------------------------------------------------------')
disp(' ')

T           =  [1 0 -2 -1;...
                0 1 2 0; ....
                0 0 1 1;...
                0 0 0 1];
Tinv        = inv(T)

%% Problem 2

%%% Social Network A-------------------------------------------------------
Aa          =  [1/4 1/4 1/4 1/4;...
                1/3 1/3 1/3 0;
                1/4 1/4 1/4 1/4;
                1/3 0 1/3 1/3];
% Eigenvalues
lambda_a    = eig(Aa);

disp('--Problem 2a ------------------------------------------------------')
disp(' ')

% Check Algebraic and Geometric Multiplcities
AlgMulta    = GetAlgebraicMultiplicity(lambda_a);
GeoMulta    = GetGeometricMultiplicity(lambda_a, Aa);
[Ta,Ja]     = CheckDiagonalizable(AlgMulta, GeoMulta, lambda_a, Aa)

% left evecs are rows of T^-1
Tainv       = inv(Ta)

%%% Social Network B-------------------------------------------------------
Ab          =  [1/3 1/3 0 1/3;...
                1/3 1/3 1/3 0;
                0 1/3 1/3 1/3;
                1/3 0 1/3 1/3];
% Eigenvalues
lambda_b    = eig(Ab);

disp('--Problem 2b ------------------------------------------------------')
disp(' ')

% Check Algebraic and Geometric Multiplcities
AlgMultb    = GetAlgebraicMultiplicity(lambda_b);
GeoMultb    = GetGeometricMultiplicity(lambda_b, Ab);
[Tb,Jb]     = CheckDiagonalizable(AlgMultb, GeoMultb, lambda_b, Ab)

% left evecs are rows of T^-1
Tbinv       = inv(Tb)

%%% Social Network C-------------------------------------------------------
Ac          =  [1/3 1/3 0 1/3;...
                1/3 1/3 1/3 0;
                0 1/2 1/2 0;
                1/2 0 0 1/2];
% Eigenvalues
lambda_c    = eig(Ac);

disp('--Problem 2c ------------------------------------------------------')
disp(' ')

% Check Algebraic and Geometric Multiplcities
AlgMultc    = GetAlgebraicMultiplicity(lambda_c);
GeoMultc    = GetGeometricMultiplicity(lambda_c, Ac);
[Tc,Jc]     = CheckDiagonalizable(AlgMultc, GeoMultc, lambda_c, Ac)

% left evecs are rows of T^-1
Tcinv       = inv(Tc)

%%% Social Network D-------------------------------------------------------
Ad          =  [1/2 0 0 1/2;...
                0 1/2 1/2 0;
                0 1/2 1/2 0;
                1/2 0 0 1/2];
% Eigenvalues
lambda_d    = eig(Ad);

disp('--Problem 2d ------------------------------------------------------')
disp(' ')

% Check Algebraic and Geometric Multiplcities
AlgMultd    = GetAlgebraicMultiplicity(lambda_d);
GeoMultd    = GetGeometricMultiplicity(lambda_d, Ad);
[Td,Jd]     = CheckDiagonalizable(AlgMultd, GeoMultd, lambda_d, Ad)

% left evecs are rows of T^-1
Tdinv       = inv(Td)

%% Problem 3
disp('--Problem 3 --------------------------------------------------------')
disp(' ')

% Eigenvectors - from hand written
V1  = [-1; 0; 1]*1/sqrt(2);
V2  = [1;-1;1]*1/sqrt(3);
V3  = [1;1;1]*1/sqrt(3);

T   = [V1,V2,V3];
Tinv= inv(T)


%%% Functions--------------------------------------------------------------

% Function to find Alg Mult of each Eigenvalue
function am = GetAlgebraicMultiplicity(evals)
    % Initialize
    lam_uni = uniquetol(evals,1e-8); % unique evals
    am      = zeros(size(lam_uni));

    % Check to see how many times each eigenvalue repeats
    for i = 1:length(am)
        ind     = find(abs(lam_uni(i) - evals) <= 1e-8);
        am(i)   = length(ind);
        fprintf('Eigenvalue of %.4f has Alg Mult of %i \n',lam_uni(i), am(i))
    end
end

% Function to find Geo Mult of each Eigenvalue
function gm = GetGeometricMultiplicity(evals, A) 
    % Initialize
    lam_uni = uniquetol(evals,1e-8); % unique evals
    gm      = zeros(size(lam_uni));

    % Check dimension of null space of lamda*I - A
    for i = 1:length(gm)
        gm(i) = size(null(lam_uni(i)*eye(size(A)) - A),2);
    end
end

function [T,J] = CheckDiagonalizable(AlgMult, GeoMult, lambda, A)
    % Check if Matrix is diagnolizable (n distinct evals or AlgMult = GeoMult)
    if (AlgMult == GeoMult)
        disp('Matrix is diagonalizable')
        [T,J] = eig(A);
    elseif ((length(uniquetol(lambda)) == length(lambda)) == 1)
        disp('Matrix is diagonalizable')
        [T,J] = eig(A);
    else
        disp('Matrix is not diagonalizable')
        [T,J] = jordan(A);
    end
end