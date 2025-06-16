
function ic = icGLM(X,Y,varargin)

% icGLM -- information criterion in Generalized Linear Models
%  Usage
%    ic = icGLM(X,Y,<optional values>)
%  Inputs
%    X            n x m real matrix; design matrix
%    Y            real vector, size n; observations (sample)
%    linkfunction char/string (default 'canonical')
%    model        char/string (default 'normal')
%    ic           char, 'AIC' (default), 'BIC'
%    submodels    cell of index vectors
%    variable name, value   
%  Outputs
%    ic            value of the information criterion
%  See also
%    help IRLS

icname = 'AIC';
[n m] = size(X);
if size(Y,1) ~= n, Y = Y'; end
if size(Y,1) ~= n, error('dimension mismatch X and Y'), end
r = size(Y,2);

% Set default values
stdev = NaN;
phi = 1;
model = 'normal';
submodels = {(1:m)};
linkfunction = 'canonical';
varname = {};
nvarargin = length(varargin);
maxit = 100;
for v=1:nvarargin,
   varv = varargin{v};
   if ischar(varv), switch varv
   case {'bernoulli','bin','binomial'},
      model = 'bin';
   case {'Poisson','poisson'}
      model = 'poisson';
   case {'normal'}
      model = 'normal';
   case {'normalstdev'},
      model = 'normalstdev';
   case {'logit'}
      linkfunction = 'logit';
   case {'probit'}
      linkfunction = 'probit';
   case {'log'}
      linkfunction = 'log';
   case {'eye','identical'},
      linkfunction = 'eye';
   case {'AIC','aic'},
      icname = 'AIC';
   case {'BIC','bic'},
      icname = 'BIC';
   case {'Cp','Mallows'},
      icname = 'Cp';
   otherwise
      varname = varv;
   end
   else
      [rv cv] = size(varv);
      if cv>1 & rv==1, varv=varv'; rv=cv; cv=1; end
      if ~isempty(varname), eval([varname '= varv;']); varname = {};
      elseif rv==1 && cv==1 && isnumeric(varv), maxit = varv;
      elseif iscell(varv) & cv==1, submodels = varv;
      end
   end
end
if isnan(stdev), stdev = phi; end
if length(linkfunction)>4, if strcmpi(linkfunction(1:5),'canon'),
switch model
case {'normal','normalstdev'},
   linkfunction = 'eye';
case {'bin'},
   linkfunction = 'logit';
case {'exp'},
   linkfunction = 'explink';
case {'Poisson','poisson','pois'}
   linkfunction = 'log';
end, end, end

ns = length(submodels);
ic = zeros(ns,1);
for s=1:ns,
  submodel = submodels{s}; p=length(submodel); Xp = X(1:n,submodel);
  betahat0 = zeros(p,1);
  [betahat muhat thetahat] = IRLS(Xp,Y,model,linkfunction,betahat0,maxit);
  logL = loglikelihood(model,thetahat,phi,Y);
  switch icname,
  case {'AIC'},
     ic(p) = 2*logL-2*p/n;
  case {'BIC'},
     ic(p) = 2*logL-log(n)*p/n;
  case {'Cp'},
     SSE = norm(muhat-Y)^2;
     ic(s) = SSE/stdev^2+2*p-n;
  end
end

% AIC = 2*logL-2*p/n;
% BIC = 2*logL-log(n)*p/n;
% Cp = SSE/stdev^2+2*p-n;

function [b d c] = bcd(model,theta,phi,y);

d = 1;
switch model,
case {'normal'},
   b = theta.^2/2; d = phi.^2; c = -y.^2/(2*phi)-log(sqrt(2*pi)*phi);
case {'bin'},
   n = 1; b = n*(theta+log(1+exp(-theta))); c=0;
   if n>1, c = log(factorial(n))-log(factorial(y))-log(factorial(n-y)); end
case {'Poisson','poisson','pois'},
   b = exp(theta); c = -log(factorial(y));
case {'exp','gamma'},
   r = 1;
   b = r*log(theta); c = (r-1)*log(y)-log(Gamma(r));
end

function logL = loglikelihood(model,theta,phi,y);

switch model,
case {'bintodo'}, % case 'bin' is now handled correctly by the general case
   n = 1; p = 1./(1+exp(-theta)); p = max(eps,min(1-eps,p));
   logL = mean(log(p).*y+log(1-p).*(1-y));
case {'normalstdev'}, % stdev is estimated within the selection procedure
   stdev2hat = mean((y-theta).^2); 
   logL = -1/2-log(2*pi*stdev2hat)/2;
otherwise,
   M = 1.e10; i = find(abs(theta)>M); theta(i) = sign(theta(i))*M;
   [b d c] = bcd(model,theta,phi,y);
   logL = mean((theta.*y-b)./d+c);
      % this expression implies the assumption that eta(theta) = theta
end
