
function [Cp AIC BIC betahats] = CpAICBIClinearmodel(Y,X,models,stdev);

% CpAICBIClinearmodel -- Computes Cp, AIC, BIC for submodels in linear model
%  Usage
%    [Cp AIC BIC betahats] = CpAICBIClinearmodel(Y,X,models,stdev);
%  Inputs
%    Y      real vector; observed response values
%    X      real rectangular matrix; full model of covariates (design matrix)
%    models cell of index vectors (integers) with unequal lengths,
%           defining submodels
%    stdev  positive real number (optional)
%  Outputs
%    Cp     real vector of length(models) (studentized version of Cp)
%    AIC    real vector of length(models)
%    BIC    real vector of length(models)
%    betahats cell of real vectors of unequal lengths
%  Note
%    This routine defines Cp = SSE/stdev^2+2*p-n;
%    For Lambda_p = SSE/n+2*p/n*stdev^2-stdev^2; multiply:
%        Lambda_p = Cp*stdev^2/n;
%  Description
%    Computes values of Cp, AIC, and BIC for all submodels defined by the
%    variable models. 
%    The full model is given by Y = X*beta + stdev*Z
%    The variance (stdev^2) is supposed to be unknown in the assessment of AIC
%    and BIC.
%  See also
%    help squareddifKLlinearmodel
%    help PElinearmodel

if nargin<4, stdev=NaN; end
n = length(Y); nmodels = length(models);
if isnan(stdev),
   I = eye(n); Pfull = X*((X'*X)\X'); pfull = size(X,2);
   stdevhatfull = norm((I-Pfull)*Y)/sqrt(n-pfull);
   stdev=stdevhatfull;
end
betahats = cell(1,nmodels);
logL = zeros(1,nmodels); SSE = zeros(1,nmodels);
p = zeros(1,nmodels);
for m = 1:nmodels,
   modelm = models{m};
   p(m) = length(modelm);
   if p(m)>0,
      Xm = X(1:n,modelm);
      betahatm = pinv(Xm)*Y;
      muhatm = Xm*betahatm;
   else
      % muhatm = mean(Y)*ones(size(Y));
      muhatm = zeros(size(Y));
      betahatm = [];
   end
   stdev2hatm = mean((Y-muhatm).^2); % MLE
   logL(m) = -1/2-log(2*pi*stdev2hatm)/2;
   SSE(m) = stdev2hatm*n;
   betahats{m} = betahatm;
   % equivalent to: SSE(m) = norm(Y-muhatm)^2;
end 
AIC = 2*logL-2*p/n;
BIC = 2*logL-log(n)*p/n;
Cp = SSE/stdev^2+2*p-n;
