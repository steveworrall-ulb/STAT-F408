
function Mp = initmom(x,p);

% initmom -- ThreshLab2/Irregular/1D/ -- 
%     initializes moments of scaling basis functions at finest scale
%  Usage
%    Mp = initmom(x,p);
%    Mp = initmom(x);
%  Inputs
%    x      grid-locations of observations (time, space)
%    p      vector of exponents (must not contain -1, default is 0)
%  Outputs
%    Mp     matrix of size length(x) times length(p)
%  Description
%    computes integral(phi(x)*x^p) with phi(x) a characteristic function on 
%    [a(i),a(i+1)], where a(1) = x(1) and a(i) = (x(i) + x(i+1))/2
%    This may serve as approximations for the finest scales scaling function
%    moments.
%  Note
%    One can use
%    fprime = initmom(f)./initmom(x); 
%    for a numerical computation of df/dx
%  See also
%    help Bsplinemoments      for moments of B-spline scaling basis functions
%    help updateprimalmoments for the construction of lifting steps on the
%                             basis of moments

if nargin < 2,
   p = 0;
end
x = row(x); n = length(x); p = column(p); nu = length(p);
Mp = zeros(nu,n);
a(n+1) = x(n); a(1) = x(1); a(2:n) = (x(1:n-1)+x(2:n))/2;
a = repmat(a,size(p));
p = repmat(p,size(x));
Mp = (a(:,2:n+1).^(p+1) - a(:,1:n).^(p+1))./(p+1);
Mp = Mp';
% alternative implementation: with integration matrix
% M = ([1 zeros(1,n-1); eye(n)]+[eye(n); zeros(1,n-1) 1])/2;
% M = M(2:n+1,1:n)-M(1:n,1:n);
% Mp = M*(a.^p)'
