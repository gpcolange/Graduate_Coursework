function F = F_NDI(q,w,I)
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

F = [(w3*(q3*w1 - q1*w3 + q4*w2))/4 - (w2*(q1*w2 - q2*w1 + q4*w3))/4 -...
    (w1*((q1*w1)/2 + (q2*w2)/2 + (q3*w3)/2))/2 +...
    (q3*(Ixz*w1^2 - Ixz*w3^2 + Ixx*w1*w3 - Izz*w1*w3))/(2*Iyy) -...
    (q2*w2*(Ixx^2*w1 + Ixz^2*w1 - Ixx*Ixz*w3 - Ixx*Iyy*w1 +...
    Ixz*Iyy*w3 - Ixz*Izz*w3))/(2*(Ixz^2 - Ixx*Izz)) +...
    (q4*w2*(Ixz^2*w3 + Izz^2*w3 - Ixx*Ixz*w1 + Ixz*Iyy*w1 -...
    Ixz*Izz*w1 - Iyy*Izz*w3))/(2*(Ixz^2 - Ixx*Izz));

    (w1*(q1*w2 - q2*w1 + q4*w3))/4 - (w3*(q2*w3 - q3*w2 + q4*w1))/4 -...
    (w2*((q1*w1)/2 + (q2*w2)/2 + (q3*w3)/2))/2 -...
    (q4*(Ixz*w1^2 - Ixz*w3^2 + Ixx*w1*w3 - Izz*w1*w3))/(2*Iyy) +...
    (q1*w2*(Ixx^2*w1 + Ixz^2*w1 - Ixx*Ixz*w3 - Ixx*Iyy*w1 +...
    Ixz*Iyy*w3 - Ixz*Izz*w3))/(2*(Ixz^2 - Ixx*Izz)) +...
    (q3*w2*(Ixz^2*w3 + Izz^2*w3 - Ixx*Ixz*w1 + Ixz*Iyy*w1 - Ixz*Izz*w1 - Iyy*Izz*w3))...
    /(2*(Ixz^2 - Ixx*Izz));

    (w2*(q2*w3 - q3*w2 + q4*w1))/4 - (w1*(q3*w1 - q1*w3 + q4*w2))/4 ...
    - (w3*((q1*w1)/2 + (q2*w2)/2 + (q3*w3)/2))/2 -...
    (q1*(Ixz*w1^2 - Ixz*w3^2 + Ixx*w1*w3 - Izz*w1*w3))/(2*Iyy) -...
    (q4*w2*(Ixx^2*w1 + Ixz^2*w1 - Ixx*Ixz*w3 - Ixx*Iyy*w1 + Ixz*Iyy*w3 -...
    Ixz*Izz*w3))/(2*(Ixz^2 - Ixx*Izz)) - (q2*w2*(Ixz^2*w3 + Izz^2*w3 -...
    Ixx*Ixz*w1 + Ixz*Iyy*w1 - Ixz*Izz*w1 - Iyy*Izz*w3))/(2*(Ixz^2 - Ixx*Izz))];