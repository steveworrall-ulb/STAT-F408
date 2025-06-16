
function [dif2 KL] = squareddifKLlinearmodel(Y,X,models,mu,stdev);

% squareddifKLlinearmodel -- Sq.diff. & KL div. for submodels in linear model
%  Usage
%    [dif2 KL] = squareddifKLlinearmodel(Y,X,models,mu,stdev);
%  Inputs
%    Y      real vector; observed response values
%    X      real rectangular matrix; full model of covariates (design matrix)
%    models cell of index vectors (integers) with unequal lengths,
%           defining submodels
%    mu     real vector; expected response values
%    stdev  positive real number (optional)
%  Outputs
%    dif2   real vector of length(models)
%    KL     real vector of length(models)
%  Description
%    Computes squared difference and Kullback-Leibler divergence for submodels
%    of a linear model Y = X*beta+stdev*Z
%    The noise in the DGP is assumed to be normal, as is the noise in all
%    submodels
%  Note
%    The expected value of dif2 is the prediction error computed in
%    PElinearmodel.
%  See also
%    help PElinearmodel
%    help CpAICBIClinearmodel


n = length(Y); nmodels = length(models);
if nargin<5, stdev = norm(Y-mu)/sqrt(n); end

ElogL = zeros(1,nmodels); dif2 = zeros(1,nmodels);
for m = 1:nmodels,
   Xm = X(1:n,models{m});
   if size(Xm,2)>0,
      betahatm = pinv(Xm)*Y;
      muhatm = Xm*betahatm;
   else
      % muhatm = mean(Y)*ones(size(Y));
      muhatm = zeros(size(Y));
   end
   stdev2hatm = mean((Y-muhatm).^2);
   ElogL(m) = -sum((mu-muhatm).^2)/(2*n*stdev2hatm)-...
                     stdev^2/(2*stdev2hatm)-log(2*pi*stdev2hatm)/2;
   dif2(m) = sum((mu-muhatm).^2)/n;
end
KL = -1/2-log(sqrt(2*pi)*stdev)-ElogL;
