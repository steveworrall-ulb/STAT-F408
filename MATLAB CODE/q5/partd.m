% Generate 10,000 samples
X = randexpinvsqrt(10000,1);

% Estimate the variance
varest = var(X);  % or: mean(X.^2) - mean(X)^2
disp(['Estimated variance: ', num2str(varest)]);

% Plot a histogram
edges = linspace(-1, 1, 41);
histogram(X, edges, 'Normalization', 'pdf'); hold on

x = linspace(-1,1,1000);        % fine grid
f = exp(-1 ./ sqrt(1 - x.^2));
dx = x(2) - x(1);
K = 1 / sum(f * dx);
fx = K * f;

plot(x, fx, 'r', 'LineWidth', 2)
legend('Samples', 'Target density (normalized)');
xlabel('x'); ylabel('Density');
title('Samples from f_X(x) = K exp(-1/sqrt(1-x^2))');
hold off