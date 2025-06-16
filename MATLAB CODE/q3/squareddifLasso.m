
function dif2 = squareddifLasso(Y,X,models,lambda,mu,stdev);

% squareddifLasso
%  Usage
%    dif2 = squareddifLasso(Y,X,models,lambda,mu,stdev);
%  Inputs
%    Y      real vector; observed response values
%    X      real rectangular matrix; full model of covariates (design matrix)
%    models cell of index vectors (integers) with unequal lengths,
%           defining submodels
%    lambda positive real vector of length(models); thresholds
%    mu     real vector; expected response values
%    stdev  positive real number (optional)
%  Outputs
%    dif2   real vector of length(models)
%  Description
%    Computes squared difference for shrinkage estimator
%    betahats = inv(Xs'*Xs)*ST(Xs'*Y) in submodels {s} of a linear model
%    Y = X*beta+stdev*Z and 
%    The noise in the DGP is assumed to be normal, as is the noise in all
%    submodels
%  See also
%    help squareddifKLlinearmodel (for orthogonal projection estimator without
%                                  shrinkage)

n = length(Y); nmodels = length(models);
if nargin<6, stdev = norm(Y-mu)/sqrt(n); end

dif2 = zeros(1,nmodels);
for m = 1:nmodels,
   Xm = X(1:n,models{m});
   if size(Xm,2)>0,
      betahatm = pinv(Xm'*Xm)*ST(Xm'*Y,lambda(m));
      muhatm = Xm*betahatm;
   else
      % muhatm = mean(Y)*ones(size(Y));
      muhatm = zeros(size(Y));
   end
   dif2(m) = sum((mu-muhatm).^2)/n;
end

