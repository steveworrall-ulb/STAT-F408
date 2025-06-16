
function void = drawgrids(xticks,yticks,varargin)

% drawgrids - ThreshLab/Plots -- 
%   Usage
%     drawgrids(xticks,yticks,plotoptions)
%   Inputs
%     xticks,yitcks   x and y values where vertical and horizontal lines will
%                     be located
%     plotoptions     (optional) plot options as in plot()
%   Outputs
%     void
%   Description
%     Adds grid of horizontal and vertical lines to a plot
%   See also
%     help drawaxes

notinlegend = {'HandleVisibility','off'};
plotoptions = {}; 
nvarargin = length(varargin); k = 1;
while k <= nvarargin, vark = varargin{k}; if ischar(vark), switch(vark)
   otherwise % the argument cannot be interpreted by this routine and it is
             % passed to the next level
      plotoptions = {plotoptions{:},vark};
      k = k+1;
   end
   else % if ischar(vark)
      plotoptions = {plotoptions{:},vark};
      k = k+1;
end, end
if length(plotoptions)<1, plotoptions = {'k-'}; end
plotoptions = {plotoptions{:},notinlegend{:}};

holdstate = ishold;
hold on
if length(yticks)>0 & length(xticks)>0,
   for h=row(yticks), plot(xticks([1 end]),[h h],plotoptions{:}), end
   for v=row(xticks), plot([v v],yticks([1 end]),plotoptions{:}), end
else
   disp('drawgrids: no grid plotted: less than 2 knots on x and/or y axis')
end
if holdstate==false, hold off, end
