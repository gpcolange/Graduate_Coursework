function u = U_NDI(q,w,I,rddot,v)
q1 = q(1);
q2 = q(2);
q3 = q(3);
q4 = q(4);
w1 = w(1);
w2 = w(2);
w3 = w(3);
Ixx = I(1,1);
Iyy = I(2,2);
Izz = I(3,3);
Ixz = -I(3,1);
r1ddot = rddot(1);
r2ddot = rddot(2);
r3ddot = rddot(3);
v1 = v(1);
v2 = v(2);
v3 = v(3);


u = [-(Ixz*q3^3*w1^2 - Ixx*q1^3*w2^2 - Ixx*q1^3*w3^2 - Ixx*q1^3*w1^2 + Ixz*q3^3*w2^2 + Ixz*q3^3*w3^2 - 4*Ixx*q1^2*r1ddot - 4*Ixx*q4^2*r1ddot + 4*Ixz*q3^2*r3ddot + 4*Ixz*q4^2*r3ddot - 4*Ixx*q1^2*v1 - 4*Ixx*q4^2*v1 + 4*Ixz*q3^2*v3 + 4*Ixz*q4^2*v3 - 4*Ixx*q1*q2*v2 - 4*Ixx*q1*q3*v3 + 4*Ixx*q2*q4*v3 - 4*Ixx*q3*q4*v2 + 4*Ixz*q1*q3*v1 - 4*Ixz*q1*q4*v2 + 4*Ixz*q2*q3*v2 + 4*Ixz*q2*q4*v1 + 2*Ixz*q4^3*w1*w2 + 2*Iyy*q4^3*w2*w3 - 2*Izz*q4^3*w2*w3 - Ixx*q1*q2^2*w1^2 - Ixx*q1*q2^2*w2^2 - Ixx*q1*q3^2*w1^2 - Ixx*q1*q2^2*w3^2 - Ixx*q1*q3^2*w2^2 - Ixx*q1*q4^2*w1^2 - Ixx*q1*q3^2*w3^2 - Ixx*q1*q4^2*w2^2 - Ixx*q1*q4^2*w3^2 + Ixz*q1^2*q3*w1^2 + Ixz*q1^2*q3*w2^2 + Ixz*q2^2*q3*w1^2 + Ixz*q1^2*q3*w3^2 + Ixz*q2^2*q3*w2^2 + Ixz*q3*q4^2*w1^2 + Ixz*q2^2*q3*w3^2 + Ixz*q3*q4^2*w2^2 + Ixz*q3*q4^2*w3^2 - 4*Ixx*q1*q2*r2ddot - 4*Ixx*q1*q3*r3ddot + 4*Ixx*q2*q4*r3ddot - 4*Ixx*q3*q4*r2ddot + 4*Ixz*q1*q3*r1ddot - 4*Ixz*q1*q4*r2ddot + 4*Ixz*q2*q3*r2ddot + 4*Ixz*q2*q4*r1ddot + 2*Ixz*q1^2*q4*w1*w2 + 2*Ixz*q2^2*q4*w1*w2 + 2*Ixz*q3^2*q4*w1*w2 + 2*Iyy*q1^2*q4*w2*w3 + 2*Iyy*q2^2*q4*w2*w3 + 2*Iyy*q3^2*q4*w2*w3 - 2*Izz*q1^2*q4*w2*w3 - 2*Izz*q2^2*q4*w2*w3 - 2*Izz*q3^2*q4*w2*w3)/(2*q4*(q1^2 + q2^2 + q3^2 + q4^2));
                                                                                                                                                                                                                                                                                                                                                                               (2*Ixz*q4^3*w1^2 - 2*Ixz*q4^3*w3^2 + Iyy*q2^3*w1^2 + Iyy*q2^3*w2^2 + Iyy*q2^3*w3^2 + 4*Iyy*q2^2*r2ddot + 4*Iyy*q4^2*r2ddot + 4*Iyy*q2^2*v2 + 4*Iyy*q4^2*v2 + 4*Iyy*q1*q2*v1 + 4*Iyy*q1*q4*v3 + 4*Iyy*q2*q3*v3 - 4*Iyy*q3*q4*v1 + 2*Ixx*q4^3*w1*w3 - 2*Izz*q4^3*w1*w3 + 2*Ixz*q1^2*q4*w1^2 + 2*Ixz*q2^2*q4*w1^2 - 2*Ixz*q1^2*q4*w3^2 + 2*Ixz*q3^2*q4*w1^2 - 2*Ixz*q2^2*q4*w3^2 - 2*Ixz*q3^2*q4*w3^2 + Iyy*q1^2*q2*w1^2 + Iyy*q1^2*q2*w2^2 + Iyy*q2*q3^2*w1^2 + Iyy*q1^2*q2*w3^2 + Iyy*q2*q3^2*w2^2 + Iyy*q2*q4^2*w1^2 + Iyy*q2*q3^2*w3^2 + Iyy*q2*q4^2*w2^2 + Iyy*q2*q4^2*w3^2 + 4*Iyy*q1*q2*r1ddot + 4*Iyy*q1*q4*r3ddot + 4*Iyy*q2*q3*r3ddot - 4*Iyy*q3*q4*r1ddot + 2*Ixx*q1^2*q4*w1*w3 + 2*Ixx*q2^2*q4*w1*w3 + 2*Ixx*q3^2*q4*w1*w3 - 2*Izz*q1^2*q4*w1*w3 - 2*Izz*q2^2*q4*w1*w3 - 2*Izz*q3^2*q4*w1*w3)/(2*q4*(q1^2 + q2^2 + q3^2 + q4^2));
 (Izz*q3^3*w1^2 - Ixz*q1^3*w2^2 - Ixz*q1^3*w3^2 - Ixz*q1^3*w1^2 + Izz*q3^3*w2^2 + Izz*q3^3*w3^2 - 4*Ixz*q1^2*r1ddot - 4*Ixz*q4^2*r1ddot + 4*Izz*q3^2*r3ddot + 4*Izz*q4^2*r3ddot - 4*Ixz*q1^2*v1 - 4*Ixz*q4^2*v1 + 4*Izz*q3^2*v3 + 4*Izz*q4^2*v3 - 4*Ixz*q1*q2*v2 - 4*Ixz*q1*q3*v3 + 4*Ixz*q2*q4*v3 - 4*Ixz*q3*q4*v2 + 4*Izz*q1*q3*v1 - 4*Izz*q1*q4*v2 + 4*Izz*q2*q3*v2 + 4*Izz*q2*q4*v1 - 2*Ixx*q4^3*w1*w2 + 2*Ixz*q4^3*w2*w3 + 2*Iyy*q4^3*w1*w2 - Ixz*q1*q2^2*w1^2 - Ixz*q1*q2^2*w2^2 - Ixz*q1*q3^2*w1^2 - Ixz*q1*q2^2*w3^2 - Ixz*q1*q3^2*w2^2 - Ixz*q1*q4^2*w1^2 - Ixz*q1*q3^2*w3^2 - Ixz*q1*q4^2*w2^2 - Ixz*q1*q4^2*w3^2 + Izz*q1^2*q3*w1^2 + Izz*q1^2*q3*w2^2 + Izz*q2^2*q3*w1^2 + Izz*q1^2*q3*w3^2 + Izz*q2^2*q3*w2^2 + Izz*q3*q4^2*w1^2 + Izz*q2^2*q3*w3^2 + Izz*q3*q4^2*w2^2 + Izz*q3*q4^2*w3^2 - 4*Ixz*q1*q2*r2ddot - 4*Ixz*q1*q3*r3ddot + 4*Ixz*q2*q4*r3ddot - 4*Ixz*q3*q4*r2ddot + 4*Izz*q1*q3*r1ddot - 4*Izz*q1*q4*r2ddot + 4*Izz*q2*q3*r2ddot + 4*Izz*q2*q4*r1ddot - 2*Ixx*q1^2*q4*w1*w2 - 2*Ixx*q2^2*q4*w1*w2 - 2*Ixx*q3^2*q4*w1*w2 + 2*Ixz*q1^2*q4*w2*w3 + 2*Ixz*q2^2*q4*w2*w3 + 2*Ixz*q3^2*q4*w2*w3 + 2*Iyy*q1^2*q4*w1*w2 + 2*Iyy*q2^2*q4*w1*w2 + 2*Iyy*q3^2*q4*w1*w2)/(2*q4*(q1^2 + q2^2 + q3^2 + q4^2))];
end