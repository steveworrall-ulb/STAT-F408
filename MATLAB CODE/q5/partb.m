n = 10000;      % Number of random values to generate
a = 2; 
b = 2;
S = -log(rand(a+b, n));           % (a+b) rows of n iid exp(1) RVs
U = sum(S(1:a, 1:n)) ./ sum(S);   % For each column: sum first a, divide by total sum
V = 2*U - 1;                      % Map from [0,1] to [-1,1]

% Plot a histogram to verify the shape matches g(x) = L*(1-x^2)
edges = linspace(-1,1,41);
histogram(V, edges, 'Normalization', 'pdf')
hold on
x = linspace(-1,1,1000);
plot(x, (3/4)*(1-x.^2), 'r','LineWidth',2)
legend('Empirical','(3/4)*(1-x^2)');
title('Random samples V and proposal density g_X(x)');
xlabel('x'); ylabel('Density');
hold off
