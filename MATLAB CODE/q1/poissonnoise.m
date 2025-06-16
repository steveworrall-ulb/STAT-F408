
function X = poissonnoise(mu,realrandom);

% poissonnoise -- generates (pseudo) Poisson count data for given intensities
%  Usage
%    X = poissonnoise(mu,realrandom);
%    X = poissonnoise(mu);
%  Inputs
%    mu     vector or matrix with positive or zero real numbers
%    realrandom binary/logical variable; default false
%  Outputs
%    X      vector or matrix with integer values; size(mu)
%  Description
%    X contains pseudo Poisson count data with intensities in mu
%  Note
%    poissonnoise calls matlab's RAND function and therefore changes RAND's
%    state.
%  See also
%    help gammanoise

if nargin < 2 | ~realrandom,
   rand('seed',931316785);
end
X = zeros(size(mu));
for i = 1:size(mu,1),
   for j = 1:size(mu,2),
     X(i,j) = 0;
     sumarrivaltime = -log(rand);
     while sumarrivaltime < mu(i,j)
       sumarrivaltime = sumarrivaltime - log (rand);
       X(i,j) = X(i,j) + 1;
     end
   end
end

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material. 
