function cross_mat      =  skew(V)
%%% Determines the cross product equivalent matrix of 3x1 vector
%%% [Vx] = [0 -V3 V2; V3 0 -V1; -V2 V1 0]

if length(V) ~= 3
    error('Input must be 3x1 Vector')
end

cross_mat   = [ 0 -V(3) V(2);...
                V(3) 0 -V(1);...
                -V(2) V(1) 0];
end