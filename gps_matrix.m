%Generates a pulse-shaping matrix with gaussian filter taps
function [G] = gps_matrix(N, N_s)
    %Generate filter taps
    g = gausswin(N_s);

    %Normalize
    g = g/sqrt(g'*g);

    %Construct matrix
    G = kron(eye(N), g);
end