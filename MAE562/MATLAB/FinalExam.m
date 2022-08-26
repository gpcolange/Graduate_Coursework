%% MAE562 Final Exam Problem 4

clear
clc

B_V         = [100 10 -10]';                                                    % Intertial velocity expressed in body frame 
B_w         = [0.1 0 0.1]';                                                     % Angular velocity of B wrt I expressed in body frame
r           = [30 20 -120]';                                                    % Position vector expressed in inertial frame

%%% 3-2-1 Euler Angle
psi         = 5;                                                                % Rotation about 3 axis [deg]
theta       = 10;                                                               % Rotation about 2 axis [deg]
phi         = 5;                                                                % Rotation about 1 axis [deg]

IcA         = [cosd(psi) -sind(psi) 0; sind(psi) cosd(psi) 0; 0 0 1];           % Rotation matrix from Inertial to first intermediate frame
AcA2        = [cosd(theta) 0 sind(theta); 0 1 0; -sind(theta) 0 cosd(theta)];   % Rotation matrix from first intermediate frame to second intermediate frame
A2cB        = [1 0 0; 0 cosd(phi) -sind(phi); 0 sind(phi) cosd(phi)];           % Rotation matrix second intermediate frame to body frame

IcB         = IcA*AcA2*A2cB;                                                    % 3-2-1 Euler Rotation Matrix

I_V         = IcB*B_V;                                                          % Intertial velocity expressed in inertial frame 
I_w         = IcB*B_w;                                                          % Angular velocity of B wrt I expressed in inertial frame

%% Spatial Velocity

skew_w_I    = [ 0 -I_w(3) I_w(2); I_w(3) 0 -I_w(1); -I_w(2) I_w(1) 0];          % Skew symmetric matrix of angular velocity of B wrt I expressed in inertial frame

I_V_B       = [skew_w_I, -cross(I_w,r)+I_V; zeros(1,4)];                        % Spatial Velocity

disp('The spatial velocity matrix is: ')
disp(I_V_B)