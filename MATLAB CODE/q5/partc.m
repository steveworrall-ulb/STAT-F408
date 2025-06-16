% Number of desired random samples
n = 10000;

% Precompute normalizations using a fine x grid
x_grid = linspace(-1,1,10000);
dx = x_grid(2) - x_grid(1);
f_grid = exp(-1 ./ sqrt(1 - x_grid.^2));
g_grid = 1 - x_grid.^2;

K = 1 / sum(f_grid * dx);      % Normalization for target
L = 1 / sum(g_grid * dx);      % Normalization for proposal, should be about 0.75

M = max((K * f_grid) ./ (L * g_grid));   % Envelope constant

% Vector for results, and counter for how many remain to be filled
results = zeros(1, n);
filled = 0;

while filled < n
    batch = n - filled;    % how many more samples we need
    
    % ----- Proposal step: generate 'batch' values from g_X(x) -----
    a = 2; b = 2;
    S = -log(rand(a+b, batch));                % draw from exponential
    U_beta = sum(S(1:a, :), 1) ./ sum(S, 1);   % Beta(2,2) random variables
    V = 2 * U_beta - 1;                        % map to [-1,1]
    
    % ----- Compute densities at each candidate -----
    f_V = exp(-1 ./ sqrt(1 - V.^2));
    g_V = 1 - V.^2;
    
    % ----- Generate uniform random values -----
    U = rand(1, batch);
    
    % ----- Acceptance rule -----
    accept = U < (K * f_V) ./ (M * L * g_V);
    
    % Store the accepted samples
    num_accepted = sum(accept);
    results(filled + 1 : filled + num_accepted) = V(accept);
    filled = filled + num_accepted;
end

% results now contains n samples from the target density

% ----- Optional: check with histogram -----
edges = linspace(-1, 1, 41);
histogram(results, edges, 'Normalization', 'pdf'); hold on
plot(x_grid, K * f_grid, 'r', 'LineWidth', 2);
legend('Accepted samples', 'Target density');
xlabel('x'); ylabel('Density');
title('Rejection sampling: Samples vs. f_X(x)');

% Output
disp(['Generated ', num2str(n), ' samples from the target density.']);
