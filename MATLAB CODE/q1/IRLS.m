
function [betahat muhat thetahat iter] = IRLS(X,Y,varargin)

% IRLS -- Iteratively Reweighted Least Squares for GLM estimation
%  Usage
%    [betahat muhat thetahat iter] = IRLS(X,Y,<optional values>)
%  Inputs
%    X            n x m real matrix; design matrix
%    Y            real vector, size n; observations (sample)
%    linkfunction char/string (default 'canonical')
%    model        char/string (default 'normal')
%    maxit        integer (max. # iterations; default Inf)
%    variable name, value   
%  Outputs
%    betahat      real vector, size m, estimated parameter values in GLM
%    muhat        real vector, size n, estimated noise-free response
%    thetahat     real vector, size n, estimated natural(canonical) parameter
%    iter         integer, effective number of iteration steps
%  Description
%    estimates betahat in the model f(y) = exp[(y*theta-b(theta))/d(phi)]*h(y)
%     

[n m] = size(X);
if size(Y,1) ~= n, Y = Y'; end
if size(Y,1) ~= n, error('dimension mismatch X and Y'), end
r = size(Y,2);

% Set default values
model = 'normal';
linkfunction = 'canonical';
varname = {};
maxit = Inf;
betahat = zeros(m,1);
tol = 1.e-8;
filled = false(1,5); % keeps track which variables have been filled in
                     % [betahat maxit]
nvarargin = length(varargin);
for v=1:nvarargin,
   varv = varargin{v};
   if ischar(varv), switch varv
   case {'bernoulli','bin','binomial'},
      model = 'bin';
   case {'Poisson','poisson'}
      model = 'poisson';
   case {'normal'}
      model = 'normal';
   case {'eye','identical'},
      linkfunction = 'eye';
   case {'logit'}
      linkfunction = 'logit';
   case {'log'}
      linkfunction = 'log';
   case {'probit'}
      linkfunction = 'probit';
   otherwise
      varname = varv;
   end
   else
      [rv cv] = size(varv);
      if cv>1 && rv==1, varv=varv'; rv=cv; cv=1; end
      if ~isempty(varname), eval([varname '= varv;']); varname = {};
      elseif rv==m && cv==1 && ~filled(1), betahat = varv; filled(1)=true;
      elseif rv==1 && cv==1 && ~filled(2), maxit = varv; filled(2)=true;
      end
   end
end
if length(linkfunction)>4, if strcmpi(linkfunction(1:5),'canon'),
switch model
case {'normal'},
   linkfunction = 'eye';
case {'bin'},
   linkfunction = 'logit';
case {'exp'},
   linkfunction = 'explink';
case {'Poisson','poisson','pois'}
   linkfunction = 'log';
end, end, end

XTX = X'*X; XTY = X'*Y;
iter = 0;
muhat0 = Y;
while iter < maxit,
   iter = iter+1; 
   thetahat = X*betahat; 
   muhat = invlink(thetahat,linkfunction);
   W = diag(varfromE(muhat,model));
   betahat = betahat+pinv(X'*W*X)*(XTY-X'*muhat);
   if max(max(abs(muhat-muhat0))) < tol, maxit = iter; end
   muhat0 = muhat;
end

function g = link(mu,linkfunction,r);

if nargin<2, linkfunction='eye'; end
switch linkfunction
case {'eye'},
    g = u;
case {'log'},
    g = log(u);
case {'logit'},
    if nargin<3, r = 1; end
    g = log(u./(r-u));
case {'probit'}
    g = invcumgauss(u);
case {'explink'}
    g = -1./u;
end

function u = invlink(t,linkfunction,r)

if nargin<2, linkfunction='eye'; end
switch linkfunction
case {'eye'},
    u = t;
case {'log'},
    u = exp(t);
case {'logit'},
    if nargin<3, r = 1; end
    u = r./(1+exp(-t));
case {'probit'}
    u = cumgauss(t);
case {'explink'}
    u = -1./t;
end

function v = varfromE(mu,model,phi);

% phi is a second parameter
if nargin<3, phi = ones(size(mu)); end

if nargin<2, model='normal'; end
switch model,
case {'normal'},
   v = phi;
case {'bin'},
   v = mu.*(1-mu).*phi;
case {'Poisson','poisson','pois'},
   v = mu;
case {'exp','gamma'},
   v = mu.^2.*phi;
end
