clear
clc
format long g

% Earth Parameters
mu_earth    = 398600.4415;
R_earth     = 6378.1363;

% Orbit parameters
a_1         = 4*R_earth;
i_1         = 28.5;
RAAN_1      = 60;
theta_1     = 0; % Ascending Node

% Circular orbit velocity
v1          = sqrt(mu_earth/a_1);

% Inertial Position
r1          = I_CDM_R(i_1, theta_1, RAAN_1)*[a_1;0;0]

a_2         = 6*R_earth;
i_2         = 45;
RAAN_2      = 60;
theta_2     = 120;

% Circular orbit velocity
v2          = sqrt(mu_earth/a_2);

% Inertial Position
r2          = I_CDM_R(i_2, theta_2, RAAN_2)*[a_2;0;0]

% angular momentum unit vector
h_hat       = cross(r1,r2)/norm(cross(r1,r2))

% Transfer orbit inclination and RAAN
i_T         = acosd(h_hat(3))
RAAN_T      = atan2d(h_hat(1),-h_hat(2))

% Unit vectors in orbit frames at initial and final position
r1_hat      = r1/norm(r1)
r2_hat      = r2/norm(r2)

theta1_hat  = cross(h_hat,r1_hat)
theta2_hat  = cross(h_hat,r2_hat)

% Latitude at each position
theta1_T    = atan2d(r1_hat(3),theta1_hat(3));
theta2_T    = atan2d(r2_hat(3),theta2_hat(3));

% Transfer angle
TA          = theta2_T - theta1_T

p           = 5*R_earth; % given semi latus rectum

% specific angular momentum 
h           = sqrt(mu_earth*p)

% f & g relations
f           = (1 - (norm(r2)/p)*(1 - cosd(TA)));
g           = norm(r2)*norm(r1)*sind(TA)/(h);

% get initial velocity along transfer arc
v1          = (r2 - f*r1)/g

% solve for sma
energy      = norm(v1)^2/2 - mu_earth/(norm(r1));
a           = -mu_earth/(2*energy)

% eccentriicty
e           = sqrt(1 - p/a)

% true anomaly at initial position along transfer arc
theta1_star = acosd((p/norm(r1) - 1)/e)

% fpa
gamma1      = sign(dot(r1,v1))*acosd(h/(norm(r1)*norm(v1)))

% f & g relations
fdot        = (dot(r1,v1)/(p*norm(r1)))*(1 - cosd(TA)) - ...
              (sqrt(mu_earth/p)*(1/norm(r1))*sind(TA));
gdot        = 1 - (norm(r1)/p)*(1 - cosd(TA));

% velocity at final position along transfer arc
v2          = fdot*r1 + gdot*v1

% true anomaly at final position along transfer arc
theta2_star = acosd((p/norm(r2) - 1)/e)

% fpa
gamma2      = sign(dot(r2,v2))*acosd(h/(norm(r2)*norm(v2)))

% Eccentric and mean anomalys at each position
E1          = 2*atan(tand(theta1_star/2)/sqrt((1+e)/(1-e)))
E2          = 2*atan(tand(theta2_star/2)/sqrt((1+e)/(1-e)))
M1          = E1 - e*sin(E1)
M2          = E2 - e*sin(E2)

% mean motion
n           = sqrt(mu_earth/a^3);

% time of flight
TOF         = (M2 - M1)/n


function DCM = I_CDM_R(i, theta, RAAN)
    C      = @(x) cosd(x);
    S      = @(x) sind(x);
    
   DCM     = [C(RAAN)*C(theta) - S(RAAN)*C(i)*S(theta), -C(RAAN)*S(theta) - S(RAAN)*C(i)*C(theta),...
                  S(RAAN)*S(i); S(RAAN)*C(theta) + C(RAAN)*C(i)*S(theta), -S(RAAN)*S(theta) + C(theta)*C(i)*C(RAAN),...
                  -C(RAAN)*S(i);S(i)*S(theta), S(i)*C(theta), C(i)];
end