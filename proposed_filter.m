function [A, alpha, hankel] = proposed_filter(d, N_s, N, G, lambda)
    %Generates our proposed filter for N_s >= 2*d+1
    g = G(1:N_s,1);
    %Hankel lifting operator of the pulse-shaping taps
    hankel = zeros(d+1, N_s-d);
    for d_k = 0:d
        hankel(d_k+1,:) = g(d+1-d_k:N_s-d_k);
    end
    one_vec = ones(d+1,1);
    %Optimal solution to problem P_3
    alpha = (hankel'*hankel+lambda*eye(N_s-d))\(hankel')*one_vec;
    
    %Equation 20b + Equation 25
    alpha_tilde = zeros(N_s,1);
    alpha_tilde(d+1:end) = alpha;
    %Equation 19
    N_t = N*N_s;
    A = zeros(N, N_t);
    for n = 1:N
        e_n = zeros(N,1);
        e_n(n) = 1;
        A(n,:) = kron(e_n, alpha_tilde);
    end
end

