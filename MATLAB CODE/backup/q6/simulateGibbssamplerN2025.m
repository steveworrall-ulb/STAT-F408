if exist('studentnumber') ~= 1, studentnumber = 0; end
if exist('dimension') ~= 1, dimension = 5; end
randn('state',studentnumber)
rand('state',studentnumber)
m = dimension;
n = 10000;
Z = randn(m,n);
A = rand(m,m);
Sigma = A*A';
Kappa = inv(Sigma);
%%%%%%%%%%%%%%%%%%%%% Transformation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = A*Z;
%%%%%%%%%%%%%%%%%%%%% Gibbs sampler %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = zeros(m,n);
B = cell(1,m);
for k=1:m,
   notk = setdiff(1:m,k);
   B{k} = Sigma(k,notk)/Sigma(notk,notk);
end
Xi = zeros(m,1);
for i=1:n,
   for k=1:m,
      notk = setdiff(1:m,k);
      Xi(k) = B{k}*Xi(notk)+Z(k,i)/sqrt(Kappa(k,k));
   end
   X(1:m,i) = Xi;
end
%%%%%%%%%%%%%%%%%%%% Matropolis-Hastings sampler %%%%%%%%%%%%%%%%%%%%%%
V = zeros(m,n);
x = V(1:m,1);
for i = 2:n,
   y = NaN;
   while isnan(y),
      eta = rand(m,1)*2-1; y = x+eta;
      alfa = exp(x'*Kappa/2*x-y'*Kappa/2*y);
      U = rand;
      if alfa<U, y = NaN; end
   end
   V(1:m,i) = y;
   x = y;
end
SigmaYhat = zeros(m,m);
SigmaXhat = zeros(m,m);
SigmaVhat = zeros(m,m);
for i=1:n,
   SigmaYhat = SigmaYhat+Y(1:m,i)*Y(1:m,i)';
   SigmaXhat = SigmaXhat+X(1:m,i)*X(1:m,i)';
   SigmaVhat = SigmaVhat+V(1:m,i)*V(1:m,i)';
end
SigmaYhat = SigmaYhat/n;
SigmaXhat = SigmaXhat/n;
SigmaVhat = SigmaVhat/n;
plot(Y(1,:),Y(2,:),'b.')
hold on
plot(X(1,:),X(2,:),'r.')
hold off
