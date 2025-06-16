x = linspace(-1,1,1000);
f = exp(-1 ./ sqrt(1 - x.^2));
g = 1 - x.^2;

dx = x(2)-x(1);
K = 1/sum(f*dx);           % numerically estimated normalization for f
L = 1/sum(g*dx);           % L = 0.75
fX = K*f;
gX = L*g;

ratio = fX ./ gX;
M = max(ratio);

figure;
plot(x, fX, 'b', 'LineWidth', 2); hold on;
plot(x, M*gX, 'r--', 'LineWidth', 2);
legend('f_X(x)', 'M \cdot g_X(x)');
xlabel('x'); ylabel('Density');
title(['Envelope plot, M = ', num2str(M)]);
hold off;

% Optionally, plot the ratio as well for completeness:
figure; plot(x, ratio, 'k', 'LineWidth', 2); 
ylabel('f_X(x)/g_X(x)'); xlabel('x');
title('Envelope ratio');

