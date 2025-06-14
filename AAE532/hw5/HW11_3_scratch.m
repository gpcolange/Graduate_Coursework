% Hyperbolic anomaly vector
HA      = linspace(H1,H2,5000);

% Initialize vectors
%v_I     = zeros(3,length(HA));
r_P     = zeros(2,length(HA));
% r_I     = v_I;
% v_I(:,1)= v1_I;
% r_I(:,1)= r1_I;

for i = 1:length(HA)
    % Calculate current true anomaly 
    theta_star  = 2*atand(sqrt((e+1)/(e-1))*tanh(HA(i)/2));
% 
%     % Change in true anomaly
%     delta_ta    = theta_star - theta_star1;
% 
%     % Calculate current radius
      r           = p/(1 + e*cosd(theta_star));
% 
%     % F and G Equations
%     f_loop      = 1 - (r/p)*(1 - cosd(delta_ta));
%     g_loop      = (r*r1)*sind(delta_ta)/sqrt(mu_sun*p);
% 
%     % New position vector
%     r_I(:,i)    = f_loop*r1_I + g_loop*v1_I;
% 
%     % F and G Equations
%     fdot        = (dot(r1_I,v1_I)/(p*r1))*(1 - cosd(delta_ta)) - (1/r1)*(sqrt(mu_sun/p))*sind(delta_ta);
%     gdot        = 1 - (r1/p)*(1 - cosd(delta_ta));
% 
%     % New Velocity Vector
%     v_I(:,i)    = fdot*r1_I + gdot*v1_I;

    % Position Vector in Perifocal Frame - Centered at Moon 
    r_P(:,i)  = [r*cosd(theta_star);r*sind(theta_star)];   
end

% figure
% plot3(r_I(1,:),r_I(2,:),r_I(3,:))
% xlabel('$\hat{x}$ [km]','Interpreter','latex')
% ylabel('$\hat{y}$ [km]','Interpreter','latex')
% zlabel('$\hat{z}$ [km]','Interpreter','latex')
% grid on
% hold on 
% plot3(0,0,0,'*y','MarkerSize',18)
% hold off
% legend('Orbit','Sun')
% title('Chronos Orbit')