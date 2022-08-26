%%% Break Condition
function [value, isterminal, direction] = BREAK(t,Z)
value      = (Z(1) <= 0);
isterminal = 1;   % Stop the integration
direction  = 0;
end
