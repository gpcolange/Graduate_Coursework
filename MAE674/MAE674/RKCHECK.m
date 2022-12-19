%% rk4 check
clear
clc

t = 0:.01:10;
z0 = [1000;0];
Z(:,1) = z0;

[t,z] = ode45(@diffeq,t,z0);

for i = 1:length(t)-1
Z(:,i+1) = RungeKutta4(@diffeq,t(i),Z(:,i),.01);
end

X = RungeKutta4_Batch(@diffeq,z0,t,.01);

figure
plot(t,Z(1,:) - X(1,:))

figure
plot(t,Z(2,:) - X(2,:))


function zdot = diffeq(t,z)
c = .05;
m = 10;
R = 6378100;
G = 6.673e-11;
M = 5.9742e24;

zdot(1,1) = z(2,1);
zdot(2,1) = -G*M/(R + z(1,1))^2 + (c/m)*z(2)^2;

end