function out = LambertsProblem(TOF_des, r1, r2, TA, mu)
% Solve Lamberts Problem via bisection

if TA > 180
    out.phi = 360 - TA;
else
    out.phi = TA;
end

% Calculate Chord
out.c           = sqrt(r1^2 + r2^2 - 2*r1*r2*cosd(out.phi));

% Semi perimeter
out.s           = (r1 + r2 + out.c)/2;

% Initial min sma
out.a_min       = out.s/2;
a_min           = out.a_min;

if TA > 180
   TOF_par  = (1/3)*sqrt(2/mu)*(out.s^(3/2) + (out.s - out.c)^(3/2));
else
   TOF_par  = (1/3)*sqrt(2/mu)*(out.s^(3/2) - (out.s - out.c)^(3/2));
end

if TOF_des <= TOF_par
    error('Parabolic Transfer Needed')
end

% Minimum Energy beta
out.beta_min    = 2*asin(sqrt((out.s-out.c)/(2*out.a_min)));

% Calculate Minimum energy TOF
out.TOF_min     = sqrt(out.a_min^3/mu)*(pi - sin(pi) - (out.beta_min - sin(out.beta_min)));

% Determine orbit type
if TA <= 180
    if TOF_des > out.TOF_min
        out.type = '1B';
    else
        out.type = '1A';
    end
else
    if TOF_des > out.TOF_min
        out.type = '2B';
    else
        out.type = '2A';
    end
end

% Choose initial max sma
a_max           = 10*out.s;

% Initialize TOF 
TOF             = 0;

% Initialize counter
count           = 0;

while (abs(TOF  - TOF_des) > 1e-6) 

    % Choose initial max sma
    a       = (a_min + a_max)/2;
    
    % Initial alpha and beta
    alpha0  = 2*asin(sqrt(out.s/(2*a)));
    beta0   = 2*asin(sqrt((out.s-out.c)/(2*a)));
    
    if strcmpi(out.type,'1A')
        alpha   = alpha0;
        beta    = beta0;
    elseif strcmpi(out.type,'2A')
        alpha   = alpha0;
        beta    = -beta0;
    elseif strcmpi(out.type,'2B')
        alpha   = 2*pi - alpha0;
        beta    = -beta0;
    elseif strcmpi(out.type,'1B')
        alpha   = 2*pi - alpha0;
        beta    = beta0;
    end

    % Lambert Equation
    TOF     = sqrt(a^3/mu)*((alpha - beta) - (sin(alpha) - sin(beta)));

    if strcmpi(out.type,'1B') || strcmpi(out.type,'2B')
        if TOF < TOF_des
            a_min = a;
        else
            a_max = a;
        end
    else
        if TOF > TOF_des
            a_min = a;
        else
            a_max = a;
        end
    end

    % Increment counter 
    count = count + 1;

    if count > 5000
        disp('Max iterations reached')
        break
    end

end

out.TOF     = TOF;
out.a       = a;
out.count   = count;
out.alpha   = alpha;
out.beta    = beta;
end