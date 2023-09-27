% Simulates the OTA system as described in our GLOBECOM 2023 paper.
% This version of the computation closely follows the notation and approach
% described in the paper.
function [MSE, MSE_mf, bias, bias_mf] = montecarlo_pedagogical(N_mc, K, M, N, N_s, d, G, sigma_z, lambda, x, verbose)
    %Number of samples to send N symbols
    N_t = N*N_s;
    %Shifting matrix to model random time delays
    E = shift_matrix(N_t);
    %Generate matched filter
    A_mf = G';
    %Generate proposed filter
    [A, a, ~] = proposed_filter(d, N_s, N, G, lambda);

    MSE = 0;
    MSE_mf = 0;
    bias = 0;
    bias_mf = 0;
    for i = 1:N_mc
        if mod(i,floor(N_mc/10)) == 0 && verbose == true
            i
        end
        %Goldenbaum sequence
        S = ones(K,M).*exp(1j*2*pi*rand(K,M));
        %Pulse shaping
        s = zeros(N_t, K);
        for k = 1:K
            s(:,k) = G*(sqrt(x(k,:))');
        end

        %Communicate M times with Goldenbaum sequence
        y_vec = zeros(M,N);
        y_vec_mf = zeros(M,N);
        for m = 1:M
            %Uniform distribution for time delays
            d_k = randi(d+1, K, 1)-1;
            %Apply time delay and random phase
            v = zeros(N_t,1);
            for k = 1:K
                v = v + S(k,m)*E^(d_k(k))*s(:,k);
            end
            %Add AWGN
            z = sigma_z/sqrt(2)*randn(size(v)) + 1j*sigma_z/sqrt(2)*randn(size(v));
            v = v+z;
            %Receive filtering
            y = A*v;
            y_mf = A_mf*v;
            y_vec(m,:) = y;
            y_vec_mf(m,:) = y_mf;
        end

        %Equation 10
        f_hat_prime = zeros(N,1);
        f_hat_prime_mf = zeros(N,1);
        for m=1:M
            f_hat_prime = f_hat_prime + (y_vec(m,:)').*(y_vec(m,:).');
            f_hat_prime_mf = f_hat_prime_mf + (y_vec_mf(m,:)').*(y_vec_mf(m,:).');
        end
        f_hat_prime = f_hat_prime/K/M;
        f_hat_prime_mf = f_hat_prime_mf/K/M;

        %Compensate for noise according to Equation 14
        f_hat = f_hat_prime-norm(A)^2*sigma_z^2/K;
        f_hat_mf = f_hat_prime_mf-norm(A_mf)^2*sigma_z^2/K;

        %The real symbol-level means
        f = sum(x)'/K;
        %MSE/bias calculation
        MSE = MSE + mean((f-f_hat).^2);
        MSE_mf = MSE_mf + mean((f-f_hat_mf).^2);
        bias = bias + mean(f-f_hat);
        bias_mf = bias_mf + mean(f-f_hat_mf);
    end
    MSE = MSE/N_mc;
    MSE_mf = MSE_mf/N_mc;
    bias = bias/N_mc;
    bias_mf = bias_mf/N_mc;
end

