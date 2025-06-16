
function [betahat, muhat, S, Cp, lambda, Sopt] = LARS(Y,X,stdev,kmax)

% LARS - Least Angle Regression for LASSO (beta-version; subject to change)
%  Usage
%    [betahat muhat S Cp lambda Sopt] = LARS(Y,X,stdev,kmax)
%  Inputs
%    Y      observations in sparsity model Y = X*beta+noise
%    X      model of covariates (rectangular matrix)
%    stdev  (estimated) standard deviation (default is 1; current version does
%           not estimate stdev)
%           if stdev==NaN, then the stop criterion is based on GCV instead of
%           Mallows's Cp
%    kmax   maximum number of nonzeros (default is samplesize n=length(Y))
%  Outputs
%    betahat estimator of beta using Sopt (see below) nonzero covariates
%    muhat   estimator of response mu = X*beta
%    S       ordered set of selected variables (indices in {1,..,length(beta)})
%    Cp      Mallows's Cp (or GCV) values of nested selections
%            The first element of Cp corresponds to a selection with "-1"
%            components, i.e., lambda = Inf; The second element of Cp
%            corresponds to a selection with 0 components, i.e., lambda on the
%            edge of introducing one variable. The third component of Cp
%            corresponds to a selection of 1 component, on the edge of
%            introducing a second one, and so on.
%    lambda  values of penalty parameter corresponding to S and Cp
%            The first element of lambda is infty (see values of Cp above)
%    Sopt    set of selected variables leading to minimum Cp value
%  Description
%  Notes
%     (1) see p.74 in Elements of statistical learning, second edition;
%     downloadable from http://www-stat.stanford.edu/~tibs/ElemStatLearn/
%     see also paper LARS, Ann.Stat.2004, pages 407-499
%     (2) LARS assumes the columns of X to be normalised.
%     Set: D = diag(sqrt(sum(X.^2))); X = X/D; 
%     (Note that this normalisation is different than that in iterativeST.m)
%  See also:
%    help iterativeST
%    help LARSmodified

Y = column(Y); n = length(Y);
m = size(X,2);

if nargin < 4, kmax = n; end
if nargin < 3, stdev = 1; end

ybar = mean(Y);
Y = Y-ybar;
muhat = zeros(size(Y));
colnormX = sqrt(sum(X.^2,1)); X = X*diag(1./colnormX);

k = 0;
% Cp = [Cp(-1) Cp(0) Cp(1)]
Cp = Inf;
Cpmin = Inf;
muhatopt = zeros(size(Y)); Sopt = []; opt = [];
lambda = Inf;
S = []; % Unlike in LARSmodified, S stands for both the currently active set
        % and for the selection history. Once selected, a component is supposed
        % to remain active at all times. 
Sprime = (1:m); % Complement of S in {1,2,...,m}
% while Cp(k+1) > Cp(k+2),
kmax = min(kmax,m-1);
while k < kmax+1,
   % At this point we have an active set with k elements, where the most
   % recently added component is active but has a value zero. We will now start
   % searching for a new component. The search proceeds in the equiangular
   % direction. On our way to the next activation, the value of the most
   % recently activated variable will go away from zero.
   muhat0 = muhat;
   if k>0,
      schatS = sign(chat(S));
      XS = X(1:n,S)*diag(schatS);
      GS = XS'*XS;
      AA = 1/sqrt(ones(1,k)*(pinv(GS)*ones(k,1)));
      wS = AA*(pinv(GS)*ones(k,1));
      uS = XS*wS;
   else
      chat = X'*Y; [Chat, j] = max(abs(chat));
      % Chat is the value of lambda when activating the first component (but
      % keeping it at zero for the moment)
      uS = X(1:n,j); uS = uS/norm(uS);
      AA = 1;
   end
   a = X'*uS;
   gamma = [(Chat-chat(Sprime))./(AA-a(Sprime)), ...
            (Chat+chat(Sprime))./(AA+a(Sprime))];
   gamma(gamma<0) = Inf;
   gamma = min(gamma,[],2);
   [gammahat, j] = min(gamma);
   muhat = muhat0 + gammahat*uS;
   newSprime = [Sprime(1:j-1), Sprime(j+1:m-k)]; 
   j = Sprime(j); S = [S, j]; Sprime = newSprime;
   % j is the (k+1)st active component, but for now, the value of betahat(j) is
   % still zero

   k = k+1; % k is number of active components
   Chatnew = Chat -gammahat*AA; % lower the threshold, so that the value of 
                                % betahat(j) becomes nonzero
   chat = X'*(Y-muhat); chatS = chat(S); Chat = max(abs(chatS));
   % Check: compare two calculations of Chat==Chatnew
   % check = Chat-Chatnew; tol = 1.e-10;
   % if abs(check)>tol, warning('inconsistency in computations of lambda'), end
   % Now store Chat as lambda(k+1) 
   lambda = [lambda, Chat];
   % The following lines are optional; for testing only
   % They check the KKT-conditions in the strictly non-zero values:
   % (XS'*XS)*betahatS = ST(XS'*Y,Chat), 
   % So, take 
   %  XS = X(1:n,S);
   %  betahatS = (XS'*XS)\(XS'*muhat);
   %  betahatSbis = (XS'*XS)\ST(XS'*Y,lambda(k+1));
   %  muhatbis = XS*betahatSbis;
   %  betahatSter = (XS'*XS)\(XS'*muhatbis);
   %  checkmu = max(abs(muhat-muhatbis));
   %  checkbeta = max(abs(betahatS-betahatSter))
   %  chatSbis = XS'*(Y-muhatbis); Chatbis = max(abs(chatSbis));
   %  checkchatSbis = chatSbis-chat(S)
   %  checkchatbis = max(abs(X'*(Y-muhatbis)-chat))
   %  checkChatbis = Chat-Chatbis
   % We see that checkmu is nonzero (so the computed muhat cannot be
   % reconstructed from Xs), but the residual of the computed muhat has the
   % same maximum (threshold) as the residual of the recomputed muhat
   if isnan(stdev), % use GCV
% norm(Y-muhat)^2
      Cpk = (norm(Y-muhat)^2/n)/(1-k/n)^2;
   else
      Cpk = norm(Y-muhat)^2/n+2*k/n*stdev^2-stdev^2;
   end
   if Cpk < Cpmin,
      Cpmin = Cpk;
      Sopt = S(1:k);
      % The newly added component (k+1) of beta is still zero for the moment
      muhatopt = muhat;
      opt = k;
   end
   Cp = [Cp, Cpk];
end
muhat = muhatopt+ybar;
XSopt = X(1:n,Sopt);

% HT:
% betahatSopt = pinv(XSopt)*Y;
% ST:
betahatSopt = pinv(XSopt)*muhatopt;
betahatSoptbis = pinv(XSopt'*XSopt)*ST(XSopt'*Y,lambda(opt+1));
betahat = zeros(m,1);
betahat(Sopt) = betahatSopt;
betahat = diag(colnormX)\betahat;

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material.
