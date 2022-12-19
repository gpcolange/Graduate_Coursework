function  y    = RungeKutta4_Batch(func,y0,t,dt,u)

if exist('u','var') == 1
    y               = zeros(length(y0),length(t));
    y(:,1)          = y0;
    
    for i = 1:length(t)-1
        k1          = dt*func(t,y(:,i),u(:,i));
        k2          = dt*func(t,y(:,i) + (k1/2),u(:,i));
        k3          = dt*func(t,y(:,i) + (k2/2),u(:,i));
        k4          = dt*func(t,y(:,i) + k3, u(:,i));
    
        y(:,i+1)    = y(:,i) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
else
        y           = zeros(length(y0),length(t));
        y(:,1)      = y0;
    
    for i = 1:length(t)-1
        k1          = dt*func(t,y(:,i));
        k2          = dt*func(t,y(:,i) + (k1/2));
        k3          = dt*func(t,y(:,i) + (k2/2));
        k4          = dt*func(t,y(:,i) + k3);
    
        y(:,i+1)    = y(:,i) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end