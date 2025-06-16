function [xx yy] = plotasblocks(x,y,varargin);

% plotasblocks - ThreshLab/Plots -- plot vectors as blocks (piecewise const.)
%   Usage
%     plotasblocks(x,y)
%     plotasblocks(y)
%     plotasblocks(x,y,option)
%     [xx yy] = plotasblocks(x,y,option)
%   Inputs
%     x      vector of abscis values
%     y      vector of coordinate values
%     'mid'  (optional) place
%     option plot options, as in matlab command plot
%     'vert' all following options (until further notice) apply to jumps
%     'hor'  all following options (until further notice) apply to hor. lines
%     'novert'  invisible jumps (NaN)
%   Outputs
%     void or [xx yy]
%     xx     vector corresponding to x, such that [xx,yy] is a piecewise
%            constant polyline
%     yy     vector corresponding to y, such that [xx,yy] is a piecewise
%            constant polyline
%            plotasblocks(x,y,option) = plot(xx,yy,option)
%   Description
%     plots vector as blocks, i.e., as piecewise constants, rather than as
%     piecewise linear (polylines), which is the matlab standard way
%   See also
%     help plotasspikes
%     help plotascircles
%     help staircase

options = {};
if nargin == 1,
   y = x; x = (1:length(y));
elseif nargin == 2,
   if ischar(y),
      options = {y};
      y = x; x = (1:length(y));
   end
else
   if ischar(y),
      options = {y,varargin{:}}; 
      y = x; x = (1:length(y));
   else
      options = varargin;
   end
end
% Parsing options
noptions = length(options);
mid = false;
vertoptions = {};
horoptions = {};
isvert = false;
novert = false;
nohor = false;
k = 0;
while k < noptions, k=k+1; optk = options{k};
   isplotopt = true;
   if ischar(optk), switch(optk)
   case {'mid'},
      mid = true; isplotopt = false;
   case {'jump','jumps','vert'},
      isvert = true; isplotopt = false;
   case {'hor'},
      isvert = false; isplotopt = false;
   case {'novert'},
      novert = true; isplotopt = false;
   case {'nohor'},
      nohor = true;
   end, end
   if isplotopt,
      if isvert,
         vertoptions = {vertoptions{:},optk};
      else
         horoptions = {horoptions{:},optk};
      end
   end
end
for k = 1:length(vertoptions), optk = vertoptions{k};
    if ischar(optk), switch(optk)
    case {'none'}
       novert = true;
    end, end
end

n = length(x);
[ny my] = size(y);
if ny~=n, y = y'; [ny my] = size(y); end
if n~=ny
   error('cannot plot vectors of unequal lengths')
end
if mid == true,
   xx = zeros(2*n,1);
   xmid = column(x(1:n-1)+x(2:n))/2;
   xx(1:2:2*n-1) = [x(1); xmid];
   xx(2:2:2*n) = [xmid; x(n)];
   yy = zeros(2*n,1);
   yy(1:2:2*n-1,1:my) = y;
   yy(2:2:2*n,1:my) = y;
else
   xx = zeros(2*n-1,1);
     xx(1:2:2*n-1) = column(x); xx(2:2:2*n-1) = column(x(2:n));
   yy = zeros(2*n-1,1);
     yy(1:2:2*n-1,1:my) = y; yy(2:2:2*n-1,1:my) = y(1:n-1,1:my);
end
holdstatus = ishold;
if novert==false & nohor==false & isempty(vertoptions),
   plot(xx,yy,horoptions{:})
else
   if nohor==false,
      xxx = NaN(3*n-2,1);
      yyy = NaN(3*n-2,1);
      xxx(1:3:3*n-2) = column(x); xxx(2:3:3*n-2) = column(x(2:n));
      yyy(1:3:3*n-2) = column(y); yyy(2:3:3*n-2) = column(y(1:n-1));
      plot(xxx,yyy,horoptions{:})
      hold on
   end
   if novert==false,
      xxx = NaN(3*n-2,1);
      yyy = NaN(3*n-2,1);
      xxx(1:3:3*n-2) = column(x); xxx(3:3:3*n-2) = column(x(2:n));
      yyy(1:3:3*n-2) = column(y); yyy(3:3:3*n-2) = column(y(1:n-1));
      plot(xxx,yyy,vertoptions{:})
   end
end
if holdstatus==false, hold off, end
if nargout == 0, clear xx yy; end

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material.
