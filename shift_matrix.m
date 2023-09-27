%Generates the shift matrix which is used to represent time shifts
function [E] = shift_matrix(N)
    I = eye(N);
    E = zeros(N);
    E(2:N,:) = I(1:N-1,:);
    E(1,N) = 1;
end