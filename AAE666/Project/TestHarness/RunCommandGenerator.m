clear
close all
clc

SetupModel
out = sim('TestCommandGenerator.slx');

extract_logsout(out)
q1c = q1c(:)';
q2c = q2c(:)';
q3c = q3c(:)';
q4c = q4c(:)';

q   = [q1c;q2c;q3c;q4c];

figure
plot(tout,vecnorm(q))

phi = zeros(length(q),1);
theta = phi;
psi = phi;

for i = 1:length(q)
    q1 = q1c(i);
    q2 = q2c(i);
    q3 = q3c(i);
    q4 = q4c(i);
    A = [(q1^2 - q2^2 - q3^2 + q4^2), 2*(q1*q2 + q3*q4), 2*(q1*q3 - q2*q4);...
           2*(q1*q2 - q3*q4), (-q1^2 + q2^2 - q3^2 + q4^2), 2*(q2*q3 + q1*q4);...
           2*(q1*q3 + q2*q4), 2*(q2*q3 - q1*q4), (-q1^2 - q2^2 + q3^2 + q4^2)];

    psi(i) = atand(A(1,2)/A(1,1));
    theta(i) = -asind(A(1,3));
    phi(i) = atand(A(2,3)/A(3,3));
end

figure
plot(tout,phi)

figure
plot(tout,theta)

figure
plot(tout,psi)