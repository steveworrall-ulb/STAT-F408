function X = randexpinvsqrt(m,n)

% generates X ~ K*exp[-1/sqrt(1-x^2))]


R = m*n;
A = 0.55;
X = zeros(1,R);
reject = (1:R);
a = 2; b = 2;
while R > 0,
   S = -log(rand(a+b,R)); U = sum(S(1:a,1:R))./sum(S); V = 2*U-1;
   X(reject) = V;
   U = rand(size(V));
   fXoverMgX =  ...
   reject =  ...
   R = length(reject);
end
X = reshape(X,m,n);
