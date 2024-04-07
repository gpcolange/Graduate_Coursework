clear
close all
clc

%% Problem 4

% Time delay [s]
h           = 0.25;

% Get non-linear equation roots
w           = fzero(@(w) tan(w*h) + w, 10);

% Verify non-linear equation is near 0.
fprintf('For a time delay of h = %.2f, tan(wh) + w = %.4f if w = %.4f \n'...
        ,h, tan(w*h) + w, w)

% Period of periodic solution
T           = 2*pi/w;

% Amplitude of periodic solution
a           = -8*cos(w*h)/pi;

fprintf(['If h = %.2f, the system has an approximate period of %.3f seconds with an '...
         'amplitude of %.3f \n'],h,T,a)

% Run Model
simout      = sim("Model.slx");

% Plot
figure
plot(simout.tout, simout.logsout{1}.Values.Data)
hold on
yline([a -a],'--r')
grid minor
xlabel('time [s]')
legend('Response','Approx. Amplitude')
ylabel('x')
title_str = ['Problem 4: h = ', num2str(h), ', a = ', num2str(a),...
            ', and T  = ',num2str(T),' [s]'];
title(title_str)


%% Problem 5

% PID Gains
Kp          = 1;
Kd          = 2;

% Sytem Transfer Function: Ki * G
G           = tf([1],[1 Kd Kp 0]);

% Root Locus plot
figure
rlocusplot(G)
[poles,Ki]   = rlocus(G,0:.001:5);

% Locate Largest Ki such that poles are in RHP
Ki_max       = max(Ki(max(real(poles)) < 0));

fprintf(['The largest Ki >=0, for which the closed loop system '...
         'is asymptotically stable about q = 0 is %.3f. \n'], Ki_max)

% Check poles of Closed Loop System 
disp('The poles of the closed loop system with the max Ki are: ')
CL_poles    = pole(feedback(Ki_max*G,1))
