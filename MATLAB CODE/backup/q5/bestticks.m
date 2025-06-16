
function ticks = bestticks(l,u,nticksmax);

% bestticks - ThreshLab/Plots -- 
%   Usage
%     ticks = bestticks(l,u,nticksmax);
%     ticks = bestticks(y,nticksmax);
%   Inputs
%     y          vector with plotted data (ordinate values)
%     l          lower bound 
%     u          upper bound
%     nticksmax  maximum number of ticks
%   Outputs
%     ticks   real vector of length smaller than nticksmax+1
%   Description
%     Finds good (best?) places of ticks on axes 
%   See also
%     help drawaxes


if nargin<3, nticksmax = NaN; end
if nargin<2, u=NaN; end
if length(l) > 1, nticksmax = u; u = max(l); l = min(l); end
if isnan(nticksmax), nticksmax = 9; end
if isnan(u), u = l; l = 0; end
d = (u-l)/nticksmax;
d = d*0.99;
% If d has only one significant figure, d1 should be that figure minus 1
% The following line finds d1 as the first significant figure of d
d1 = floor(d/10^floor(log10(d)));
dd = [0.2 0.5 0.5 0.5 1 1 1 1 1];
d1 = dd(d1);
d = d1*10^ceil(log10(d));
l = ceil(l/d)*d;
ticks = (l:d:u);
