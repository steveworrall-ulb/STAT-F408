
function PE = PElinearmodel(mu,X,models,stdev);

% PElinearmodel -- Computes the prediction error for submodels in linear model
%  Usage
%    PE = PElinearmodel(mu,X,models,stdev);
%  Inputs
%    mu     real vector; expected response values
%    X      real rectangular matrix; full model of covariates (design matrix)
%    models cell of index vectors (integers) with unequal lengths,
%           defining submodels
%    stdev  positive real number (default 1)
%  Outputs
%    PE     real vector of length(models)
%  Description
%    Computes the prediction error (expected squared difference) for submodels
%    of a linear model Y = X*beta+stdev*Z
%    The noise in the DGP is assumed to be normal, as is the noise in all
%    submodels
%  Note
%    The PE computed by this routine is the expected value of the squared
%    difference, computed by squareddifKLlinearmodel.
%  See also
%    help squareddifKLlinearmodel
%    help CpAICBIClinearmodel


n = length(mu); I = eye(n);
if nargin<4, stdev = 1; end

nmodels = length(models);
p = zeros(1,nmodels);
bias2 = zeros(1,nmodels);
for m = 1:nmodels,
   p(m) = length(models{m});
   Xm = X(1:n,models{m});
   if p(m)>0,
      P = Xm*pinv(Xm);
   else
      P = zeros(n,n);
   end
   bias2(m) = mu'*(I-P)*mu/n;
end
PE = bias2 + p*stdev^2/n;
