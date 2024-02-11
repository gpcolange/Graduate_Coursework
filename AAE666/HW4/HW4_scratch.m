clear
close all
clc

%% Problem 1
syms x1 x2 lambda real
xdot = [x2; -x1 + x1^3 - x2];
V  = x2^2/2 + x1^2/2 - x1^4/4 + lambda*x2*x1 + lambda*x1^2/2;
subs(jacobian(V,[x1;x2]),[x1;x2],zeros(2,1))
subs(hessian(V,[x1;x2]),[x1;x2],zeros(2,1))
Vdot = expand(jacobian(V,[x1;x2])*xdot)

% lambda = .001:.001:.999;
% det = - lambda.^2 + lambda + 1;
% 
% figure
% plot(lambda,det)

% exp = matlabFunction(subs(Vdot,x1,0.99));
% x2 = -10:.1:10;
% lambda = 0.9;
% data = exp(lambda,x2);
% 
% figure
% plot(x2,data)

%% Problem 4
% syms x1 x2 k1 k2 lambda real
% 
% xdot  = [x2; (x1 - x1^3 -k1*x1 - k2*x2)];
% 
% V = x2^2/2 + lambda*k2*x2*x1 + lambda*k2^2*x1^2/2 + k1*x1^2/2 + x1^4/4 - x1^2/2;
% P   = 1/2*[lambda*k2^2, lambda*k2; lambda*k2, 1];
% x = [x1;x2];
% simplify(V - (x'*P*x + k1*x1^2/2 + x1^4/4 - x1^2/2))
% 
% DV0 = subs(jacobian(V,[x1;x2]),[x1;x2],zeros(2,1))
% V0 = subs(V,[x1;x2],zeros(2,1))
% %D2V = subs(hessian(V,[x1;x2]),[k1;k2;lambda],[1.2;0.5;0.5])
% D2V = subs(hessian(V,[x1;x2]),k1,1 - lambda*k2^2 + lambda^2*k2^2)

% clear
% close all
% clc
% 
% % Control gains 
% k1      = 1.2;   % From Lyapunov analysis k1 > 1
% k2      = 0.5;    % From Lyapunov analysis k2 > 0
% 
% % Initial conditions 
% IC      = zeros(2,1) + 2*randn(2,10);
% 
% % sim time
% time    = (0:.005:50)';
% 
% % ODE45 solver options
% options = odeset('AbsTol',1e-8,'RelTol',1e-8);
% 
% % Loop through all IC's
% for i = 1:length(IC)
% 
%     % Closed loop system
%     [~, X_cl]     = ode45(@(t,x) ControlledDuffingSystem(t,x,[k1 k2]),...
%                     time, IC(:,i), options);
% 
%     % Generate Plots
%     title_str     = sprintf(['Duffing System with IC x_1 = %.2f & ',...
%                             ' and x_2 = %.2f \n'],IC(1,i),IC(2,i));
% 
%     figure(i)
%     subplot(211)
%     plot(time, X_cl(:,1))
%     ylabel('x_1')
%     grid minor
%     title(title_str)
%     subplot(212)
%     plot(time, X_cl(:,2))
%     ylabel('x_2')
%     grid minor
%     xlabel('Time [s]')
% 
% end
% 
% function xdot = ControlledDuffingSystem(t,x,K)
%     % Control law, u = -k1*x1 - k2*x2
%     u           = -K*x;
% 
%     % State space model
%     xdot(1,1)   = x(2,1);
%     xdot(2,1)   = x(1,1) - x(1,1)^3 + u;
% end
