clear
clc

C = [1 0 0 0; 0 0 0 1];
A = [0 0 1 0;1 0 2 0;0 1 3 0;0 0 -21 5];
n = length(A);

co = ctrb(A',C')

[r,p] = rref(co);

L_full     = [C(1,:)',A'*C(1,:)',A'^2*C(1,:)',A'^3*C(1,:)',C(2,:)',A'*C(2,:)',A'^2*C(2,:)',A'^3*C(2,:)'];
L          = zeros(size(A));
count = 1;

for i = 1:length(L_full)
    for j  = 1:length(p)
        if co(:,p(j)) == L_full(:,i)
            L(:,count) = L_full(:,i);
            count = count + 1;
        end
    end
end



%L   = L_full(:,p)

inv_L = inv(L);

q1 = inv_L(3,:);
q2 = inv_L(4,:);

T = [q1;q1*A';q1*A'^2;q2]

A = (T*A'*inv(T))';
C = C*T';
