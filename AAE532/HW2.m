clear
close all
clc

%% Problem 1a
disp('-----Start Problem 1 ------------------')
% Earth -> Moon -> s/c -> sun -> Jupiter (Collinear)

% Assume 1-D position vectors along y axis. Positive y is right 
r_sun_jupiter   = 778279959;                    % [km]
r_sun_earth     = -149597898;                   % [km]
r_earth_moon    = 384400;                       % [km]
r_sun_moon      = r_sun_earth + r_earth_moon;
r_moon_sc       = 77500;                        % [km]
r_sun_sc        = r_sun_moon + r_moon_sc; 

% Gravitational constant
G               = 6.6743e-11*(1/1000)^3;        % [km^3/kg-s^2]

% Masses
m_sc            = 130;                          % [kg]
m_earth         = 398600.4415/G;                % [kg]
m_jupiter       = 126712767.8578/G;   
m_moon          = 4902.8005821478/G;   
m_sun           = 132712440017.99/G;

% Center of Mass Calculation [km]
disp('Center of mass location')
COM             = (r_sun_jupiter*m_jupiter + r_sun_earth*m_earth + ...
                   r_sun_moon*m_moon + r_sun_sc*m_sc)/ + (...
                   m_jupiter + m_earth + m_moon + m_sc + m_sun)

%% Problem 1b

% Relative position vectors
r_earth_sc       = r_sun_sc - r_sun_earth;
r_jupiter_sc     = r_sun_sc - r_sun_jupiter;

% Accelerations on spacecraft due to bodies
disp('Acceleration on s/c due to Earth')
a_earth_sc       = -G*m_earth*r_earth_sc/(norm(r_earth_sc)^3)
disp('Acceleration on s/c due to moon')
a_moon_sc        = -G*m_moon*r_moon_sc/(norm(r_moon_sc)^3)
disp('Acceleration on s/c due to Jupiter')
a_jupiter_sc     = -G*m_jupiter*r_jupiter_sc/(norm(r_jupiter_sc)^3) 
disp('Acceleration on s/c due to Sun')
a_sun_sc         = -G*m_sun*r_sun_sc/(norm(r_sun_sc)^3) 
disp('Net Acceleration on s/c')
a_net            = a_earth_sc + a_sun_sc + a_moon_sc + a_jupiter_sc

disp('-----End Problem 1 ------------------')

%% Problem 2a
disp('-----Start Problem 2 ------------------')

% More position vectors
r_sc_sun         = -r_sun_sc;
r_moon_sun       = -r_sun_moon;
r_sc_earth       = -r_earth_sc;
r_sc_jupiter     = r_sun_jupiter - r_sun_sc;
r_moon_jupiter   = r_sun_jupiter - r_sun_moon;
r_moon_earth     = -r_earth_moon;

disp('Dominant acceleration: ')
a_dominant       = G*(m_sc + m_moon)*r_moon_sc/(norm(r_moon_sc)^3)

disp('Total direct acceleration: ')
a_direct         = G*m_sun*(r_sc_sun/norm(r_sc_sun)^3) + ...
                   G*m_jupiter*(r_sc_jupiter/norm(r_sc_jupiter)^3) + ...
                   G*m_earth*(r_sc_earth/norm(r_sc_earth)^3)

disp('Total indirect acceleration: ')
a_indirect       = G*m_sun*(r_moon_sun/norm(r_moon_sun)^3) + ...
                   G*m_jupiter*(r_moon_jupiter/norm(r_moon_jupiter)^3) + ...
                   G*m_earth*(r_moon_earth/norm(r_moon_earth)^3)

disp('Total pertubation acceleration: ')
a_total_pert     = a_direct - a_indirect

disp('Pertubation acceleration due to Earth: ')
a_pert_earth     = G*m_earth*(r_sc_earth/norm(r_sc_earth)^3) - G*m_earth*(r_moon_earth/norm(r_moon_earth)^3)
disp('Pertubation acceleration due to Jupiter: ')
a_pert_jupiter   = G*m_jupiter*(r_sc_jupiter/norm(r_sc_jupiter)^3) -  G*m_jupiter*(r_moon_jupiter/norm(r_moon_jupiter)^3)
disp('Pertubation acceleration due to Sun: ')
a_pert_sun       = G*m_sun*(r_sc_sun/norm(r_sc_sun)^3) - G*m_sun*(r_moon_sun/norm(r_moon_sun)^3)
disp('-----End Problem 2 ------------------')

%% Problem 3a
disp('-----Start Problem 3 ------------------')

Gm_charon       = 119.480; % [km^3/s^2]
Gm_pluto        = 981.601; % [km^3/s^2]
r_pc            = 19596;   % Positive xhat direction [km]

disp('Pluto, Charon Center of Mass: ')
COM_pc          = Gm_charon*r_pc/(Gm_charon + Gm_pluto)

disp('Position vector from Center of Mass to Pluto: ')
r_cm_pluto      = -COM_pc % Negative xhat direction
disp('Position vector from Center of Mass to Charon: ')
r_cm_charon     = r_pc + r_cm_pluto

%% Problem 3b 

% Inertial Velocities along y axis [km/s]
rdot_charon     = .211319; 
rdot_pluto      = -.025717;

% Inertial velocity in y axis
rdot_pc         = rdot_charon - rdot_pluto; % [km/s]

% Specific Angular Momentum
h               = cross([r_pc;0;0],[0;rdot_pc;0]);

% Angular Velocity
disp('Angular Velocity of Charon relative to Pluto: ')
thetadot        = norm(h)/norm(r_pc)^2

%% Problem 3c

disp('Linear momentum of the system: ')
p               = Gm_pluto*rdot_pluto/G + Gm_charon*rdot_charon/G

disp('Velocity of center of mass: ')
v_cm            = p/(Gm_pluto/G + Gm_charon/G)

%% Problem 3d 
disp('C3 is: ')
C3              = Gm_pluto*r_cm_pluto*rdot_pluto/G + Gm_charon*r_cm_charon*rdot_charon/G

%% Problem 3e
T               = 1/2*(Gm_pluto/G)*rdot_pluto^2 + 1/2*(Gm_charon/G)*rdot_charon^2;
U               = G*(Gm_pluto/G)*(Gm_charon/G)/r_pc;

disp('C4 is: ')
C4              = T - U