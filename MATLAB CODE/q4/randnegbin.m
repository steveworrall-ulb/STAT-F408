
function ans = randnegbin(m,n,p,alfa)

% randnegbin -- generates (pseudo) random numbers according to the negative
%               binomial distribution (number of failures before alfath success)
%  Usage
%    ans = randnegbin(m,n,p,alfa)
%    ans = randnegbin(dim,p,alfa)
%  Inputs
%    dim    vector of integers with size of output vector: size(ans) = dim
%    m n    two integeres such that size(ans) = [m n]
%    p      real in [0,1], default = 1/2;
%           parameter of the negative binomial distribution; default = 1/2
%           probability of success in Bernoulli experiment
%    alfa   positive real
%           parameter of the negative binomial distribution; default = 1
%           (with alfa = 1, the negative binomial is the geometrical distr.)
%           when alfa is integer, this parameter is interpreted as the aimed
%           number of successes in a Bernoulli experiment
%  Outputs
%    ans    random vector of size dim or m times n
%  Description
%    generates real random numbers according to the discrete distribution
%    pdf: p(x) = Gamma(x+alfa)*(1-p)^x/x!*p^alfa/Gamma(alfa)
%    If alfa is integer, this models the number of failures before alfa 
%    successes in a Bernoulli experiment
%    An alternative definition of the negative binomial distribution takes the
%    number of trials (successes+failures) in a Bernoulli experiment as random 
%    variable; to generate values from that distribution, use
%    randnegbin(m,n,p,alfa)+alfa
%  Note
%    randnegbin calls matlab's RAND function and therefore changes RAND's
%    state.
%  Examples
%    X = randnegbin(2,4,0.3,5.3)
%    X = randnegbin(size(a),0.3,5.3)
%  See also
%    help pmfnegbin

if (nargin == 0)
  dim = 1; p = 0.5; alfa = 1;
elseif (nargin == 1)
  if (length(m) > 1)
     dim = m; p = 0.5; alfa = 1;
  else
     dim = 1; alfa = 1; p = m;
  end
elseif (nargin == 2)
  if (length(m) > 1)
     dim = m; p = n; alfa = 1;
  else
     dim = [m n]; p = 0.5; alfa = 1;
  end
elseif (nargin == 3)
  if (length(m) > 1)
     dim = m; alfa = p; p = n;
  else
     dim = [m n]; alfa = 1;
  end
else
  dim = [m n];
end

m = dim(1); n = dim(2);
mn = m*n;

lambda = log(1./(1-p));
r = ceil(alfa);
ralfa = r-alfa;
X = ceil(randexp(r,mn,lambda));
N = sum(X)-r;

if ralfa > eps,

   k = (1:mn);
   K = mn;

   while K > 0, % rejection sampling
      U = rand(1,K);
      L = pmfnegbin(N(k),p,alfa)./pmfnegbin(N(k),p,r)*p^ralfa;
      reject = find(U>L);
      k = k(reject);
      K = length(k);
      X = ceil(randexp(r,K,lambda));
      N(k) = sum(X)-r;
   end
end
ans = reshape(N,m,n);

% Copyright (c) Maarten Jansen
% 
% This software is part of ThreshLab and is copyrighted material. 
