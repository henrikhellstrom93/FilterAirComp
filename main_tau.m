clear
clc

%Number of Monte-Carlo trials 
N_mc = 10000;
%If verbose == true: Print progress regularly
verbose = true;
%Number of user devices
K = 100;
%Length of Goldenbaum sequence
M = 10;
%Number of symbols
N = 10;
%Transmitted messages (N messages per device)
x = 3*rand(K,N);
%Standard deviation of AWGN
sigma_z = 1;
%Regularization parameter. Higher lambda => more noise suppression
lambda = 0.1;
%List of maximum delays to consider for numerical results
d_list = 0:10;
%Number of samples per symbol = 2*d + N_e. Must be >= 1
N_e = 2;
%MSE of proposed method
MSE = zeros(size(d_list));
%MSE of matched filter
MSE_mf = zeros(size(d_list));
%Bias of proposed method
bias = zeros(size(d_list));
%Bias of matched filter
bias_mf = zeros(size(d_list));

tic
for d = d_list
    d
    N_s = 2*d+N_e;
    %Generate pulse-shaping filter
    G = gps_matrix(N,N_s);
    %Run m_c Monte-Carlo trials
    %[MSE_tau, MSE_mf_tau, bias_tau, bias_mf_tau] = montecarlo_pedagogical(N_mc, K, M, N, N_s, d, G, sigma_z, lambda, x, verbose);
    [MSE_tau, MSE_mf_tau, bias_tau, bias_mf_tau] = montecarlo_efficient(N_mc, K, M, N, N_s, d, G, sigma_z, lambda, x, verbose);
    %Store results
    MSE(d+1) = MSE_tau;
    MSE_mf(d+1) = MSE_mf_tau;
    bias(d+1) = bias_tau;
    bias_mf(d+1) = bias_mf_tau;
end
toc
disp("Done!")

save("random")

%% ---- Plotting ---- %%
variance = false; %Plot variance or MSE

if variance == false
    plot(d_list, MSE, 'r')
    hold on;
    plot(d_list, MSE_mf, "g")
    hold on;
else
    plot(d_list, MSE-bias.^2, 'r')
    hold on;
    plot(d_list, MSE_mf-bias_mf.^2, "g")
    hold on;
end
p = plot(d_list, bias.^2);
p.Color = "#A2142F";
hold on;
p2 = plot(d_list, bias_mf.^2);
p2.Color = "#77AC30";
if variance == false
    legend("MSE", "MF MSE", "bias^2", "MF bias^2")
    ylabel("MSE/bias^2")
else
    legend("Variance", "MF Variance", "bias^2", "MF bias^2")
    ylabel("variance/bias^2")
end
xlabel("d")
ylim([0 max(MSE_mf)+0.1])
grid on;
