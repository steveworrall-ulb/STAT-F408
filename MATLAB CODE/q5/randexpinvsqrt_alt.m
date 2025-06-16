function X = randexpinvsqrt(m,n)
%RANDEXPINVSQRT Generate random numbers with density proportional to exp(-1/sqrt(1-x^2)) on [-1,1]
% Usage:
%   X = randexpinvsqrt(m,n)
% Returns an [m x n] matrix of samples from the target density

    % Precompute normalization and envelope using a fine grid
    x_grid = linspace(-1,1,10000);
    dx = x_grid(2) - x_grid(1);
    f_grid = exp(-1 ./ sqrt(1 - x_grid.^2));
    g_grid = 1 - x_grid.^2;
    K = 1 / sum(f_grid * dx);                  % target norm
    L = 1 / sum(g_grid * dx);                  % proposal norm
    M = max((K * f_grid) ./ (L * g_grid));     % envelope

    n_samples = m*n;
    X = zeros(1, n_samples);
    filled = 0;

    a = 2; b = 2; % for proposal

    while filled < n_samples
        batch = n_samples - filled;

        % Step 1: Generate proposals from g_X(x) (Beta(2,2) mapped to [-1,1])
        S = -log(rand(a+b, batch));
        U_beta = sum(S(1:a,:), 1) ./ sum(S,1);
        V = 2*U_beta - 1;

        % Step 2: Evaluate densities
        fV = exp(-1 ./ sqrt(1 - V.^2));
        gV = 1 - V.^2;

        % Step 3: Generate uniforms for acceptance
        U = rand(1, batch);

        % Step 4: Acceptance criterion
        accept = U < (K * fV) ./ (M * L * gV);

        % Step 5: Store accepted values
        num_accepted = sum(accept);
        X(filled + 1 : filled + num_accepted) = V(accept);
        filled = filled + num_accepted;
    end

    X = reshape(X, m, n);

end
