clear
clc

%% Feasability 
% A               = randn(3,3);
% setlmis([]);
% P               = lmivar(1,[3,1]);
% lmiterm([1 1 1 P],A',1,'s');
% lmis            = getlmis;
% 
% [tmin, xfeas]   = feasp(lmis);
% P               = dec2mat(lmis, xfeas, P)

%% State Feedback Controller 
% 
% A = [0 1 0 0;0 0 -1 0;0 0 0 1; 0 0 11 0];
% B = [0;.1;0;-.1];
% [n, m] = size(B);
% 
% 
% cvx_begin sdp 
% 
% % Variable definition
% variable S(n, n) symmetric
% variable Z(m, n)
% 
% % LMIs
% S*A'-Z'*B'+A*S-B*Z+2*S <= -eps*eye(n);
% S >= eps*eye(n);
% cvx_end
% 
% K = Z/S; % compute K matrix


%% CARE
A = sqrt(2);
B = 3;
Q = 1/4;
R = 9;

[K,P] = lqr(A,B,Q,R)

A'*P + P*A + Q - P*B*inv(R)*B'*P