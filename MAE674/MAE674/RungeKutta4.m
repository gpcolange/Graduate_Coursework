function [y_kp1] = RungeKutta4(diffeq,tk,yk,dt,u)
% Performs Runge Kutta 4th order numerical intergration at time t,
% primarily to be used for EKF implementation, as it is a single step
% evaluation. y_(k+1) = RungeKutta4(y_k,dt). Get y_k by evaluating function
% at current time step

if ~exist('u','var')
    k1      = dt.*diffeq(tk,yk);
    k2      = dt.*diffeq(tk +.5, yk + (k1*.5));
    k3      = dt.*diffeq(tk +.5, yk + (k2*.5));
    k4      = dt.*diffeq(tk + 1, yk + k3);

    y_kp1   = yk + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
else
    k1      = dt.*diffeq(tk,yk,u);
    k2      = dt.*diffeq(tk +.5, yk + (k1*.5),u);
    k3      = dt.*diffeq(tk +.5, yk + (k2*.5),u);
    k4      = dt.*diffeq(tk + 1, yk + k3,u);

    y_kp1   = yk + (1/6)*(k1 + 2*k2 + 2*k3 + k4);

end

end