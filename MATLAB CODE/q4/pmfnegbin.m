function pN = pmfnegbin(n,p,alfa)

% pmfnegbin -- probability mass function of the negative binomial distribution
%              (number of failures before alfath success)
%  Usage
%    pN = pmfnegbin(n,p,alfa);
%  Inputs
%    n      matrix of integers
%    p      real in [0,1], default = 1/2;
%           parameter of the negative binomial distribution; default = 1/2
%           probability of success in Bernoulli experiment
%    alfa   positive real
%           parameter of the negative binomial distribution; default = 1
%           (with alfa = 1, the negative binomial is the geometrical distr.)
%           when alfa is integer, this parameter is interpreted as the aimed
%           number of successes in a Bernoulli experiment
%  Outputs
%    pN     matrix of size(n)
%  Description
%    Computes P(N=n), where N is number of failures before alfa successes in
%    Bernoulli experiment (if alfa is integer; generalisation otherwise)
%    If X is the number trials (failures and successes) before alfa successes,
%    then P(X=x) = P(N+alfa=x) = P(N=x-alfa)

pN = (1-p).^n./factorial(n).*gamma(n+alfa)*p^alfa/gamma(alfa);
