
function void = drawaxes(varargin);

% drawaxes - ThreshLab/Plots -- 
%   Usage
%     drawaxes
%     drawaxes(...various inputs...)
%   Inputs
%     void or optional inputs, including:
%     'color',color             character or [r g b] vector (r,g,b in [0,1])
%     'linewidth',linewidth     integer number (default is 3)
%     'arrowlength',arrowlength relative length of arrow (default is 0.07)
%     'arrowwidth',arrowwidth   relative width of arrow (default is 0.02)
%     'xtick',xtick             vector of 'tick' locations on x-axis
%     'xticklabel',xticklabel   
%     'ytick',ytick             vector of 'tick' locations on y-axis
%     'yticklabel',yticklabel   
%   Outputs
%     void
%   Description
%     Draws arrowed axes on current figure
%   See also
%     help drawgrids

notinlegend = {'HandleVisibility','off'};
% First, tag the axes
fignum = gcf;
if ~isnumeric(fignum), fignum = get(gcf,'number'); end
tag = ['drawaxes' num2str(fignum)];
% now delete existing axes with the same tag
% findobj will not find objects with {'HandleVisibility','off'}
delete(findall(0,'tag',tag))

nvarargin = length(varargin);
linecolor = 'k';
linewidth = 3;
arrowlength = 0.07;
arrowwidth = 0.02;
% plotoptions = {};
plotoptions = notinlegend;
xticks = [];
xticklabels = NaN;
yticks = [];
yticklabels = NaN;
fontsize = 12;
fontweight = 'normal';
textcolor = NaN;
interpreter = [];

k = 1;
while k <= nvarargin, vark = varargin{k}; if ischar(vark), switch(vark)
   case {'color','linecolor'},
      linecolor = varargin{k+1};
      k = k+2;
   case {'b','g','r','c','m','y','k','w'},
      linecolor = vark;
      k = k+1;
   case {'linewidth'},
      linewidth = varargin{k+1};
      k = k+2;
   case {'arrowlength'},
      arrowlength = varargin{k+1};
      k = k+2;
   case {'arrowwidth'},
      arrowwidth  = varargin{k+1};
      k = k+2;
   case {'xtick','XTick','xticks','Xtick','Xticks','XTicks'},
      xticks = varargin{k+1};
      k = k+2;
   case {'xticklabel','XTickLabel','xticklabels'},
      xticklabels = varargin{k+1};
      k = k+2;
   case {'ytick','YTick','yticks','Ytick','Yticks','YTicks'},
      yticks = varargin{k+1};
      k = k+2;
   case {'yticklabel','YTickLabel','yticklabels'},
      yticklabels = varargin{k+1};
      k = k+2;
   case {'fontsize','textsize','FontSize','TextSize','Fontsize'},
      fontsize = varargin{k+1};
      k = k+2;
   case {'fontweight','FontWeight','Fontweight'},
      fontweight = varargin{k+1};
      k = k+2;
   case {'textcolor','fontcolor'},
      textcolor = varargin{k+1};
      k = k+2;
   case {'interpreter'},
      interpreter = varargin{k+1};
      k = k+2;
   otherwise % the argument cannot be interpreted by this routine and it is
             % passed to the next level
      plotoptions = {plotoptions{:},vark};
      k = k+1;
   end
   else % if ischar(vark)
      plotoptions = {plotoptions{:},vark};
      k = k+1;
end, end
plotoptions = {plotoptions{:},'linewidth',linewidth,'color',linecolor};
if isnan(textcolor), textcolor = linecolor; end
   
figpos = get(gcf,'position'); figwidth = figpos(3); figheight = figpos(4);
Wx = arrowwidth*420/figheight;
Wy = arrowwidth*560/figwidth;
Lx = arrowlength*560/figwidth;
Ly = arrowlength*420/figheight;
holdon = ishold;
axy = axis;
x0 = axy(1); x1 = axy(2); y0 = axy(3); y1 = axy(4);
if x0*x1 > 0, Ypos = x0; else, Ypos = 0; end
if y0*y1 > 0, Xpos = y0; else, Xpos = 0; end
% Xpos: position (height) of x-axis
% Ypos: position (x-coordinate) of y-axis

Dx = x1-x0; Dy = y1-y0;
marginfactor = 1; % for use in octave;
if isoctave, marginfactor = 1.1; end
                  % without a margin, octave's axis ruins the arrows
axy(1) = min(axy(1),Ypos-marginfactor*Dx*Wy);
axy(2) = max(axy(2),Ypos+marginfactor*Dx*Wy);
axy(3) = min(axy(3),Xpos-marginfactor*Dy*Wx);
axy(4) = max(axy(4),Xpos+marginfactor*Dy*Wx);


% First, the axes
hold on
h = plot(axy(1:2)-[0 1]*Dx*Lx,Xpos*[1 1],plotoptions{:});
set(h,'tag',tag)
h = plot(Ypos*[1 1],axy(3:4)-[0 1]*Dy*Ly,plotoptions{:});
set(h,'tag',tag)

% Second, the arrows
Dx = axy(2)-axy(1);
Dy = axy(4)-axy(3);
h = fill(axy(2)-Dx*Lx*[1 0 1],Xpos+Wx*Dy*[-1 0 1],linecolor,notinlegend{:});
set(h,'tag',tag)
h = fill(Ypos+Wy*Dx*[-1 0 1],axy(4)-Dy*Ly*[1 0 1],linecolor,notinlegend{:});
set(h,'tag',tag)

% Third, the ticks
fontoptions = {'fontsize',fontsize,'fontweight',fontweight,'color',textcolor};
if ~isempty(interpreter),
   fontoptions = {fontoptions{:},'interpreter',interpreter};
end
if isnumeric(xticklabels), if all(isnan(xticklabels)), 
   xticklabels = xticks;
end, end
Nxlabels = length(xticklabels);
if isnumeric(yticklabels), if all(isnan(yticklabels)),
   yticklabels = yticks;
end, end
Nylabels = length(yticklabels);
for l = 1:length(xticks),
   xt = xticks(l);
   if Nxlabels < l,
      xl = '';
   elseif iscell(xticklabels),
      xl = xticklabels{l};
   else
      xl = xticklabels(l);
   end
   if xt < axy(2)-Dx*Lx & xt > axy(1),
      if isnumeric(xl),
         if isnan(xl), xl = ''; else, xl = num2str(xl); end
      end
      xll = length(xl);
      h = plot([xt xt],[Xpos Xpos+8/figheight*Dy],plotoptions{:});
      set(h,'tag',tag)
      lcr = 'center';
      if abs(xt-Xpos) < 0.01*Dx, xt = Xpos+0.005*Dx; lcr = 'left'; end
      h = text(xt,Xpos,xl,...
           'HorizontalAlignment',lcr,'Verticalalignment','top',fontoptions{:});
      set(h,'tag',tag)
   end
end

useformat = false;
if all(isnumeric(yticklabels)),
   if length(yticklabels) > 1,
      d = ceil(-log10(yticklabels(2)-yticklabels(1)));
      % d is required precision
      if d>0, useformat = true; ylformat = ['%.' num2str(d) 'f']; end
   end
end
for l = 1:length(yticks),
   yt = yticks(l);
   if Nylabels < l,
      yl = '';
   elseif iscell(yticklabels),
      yl = yticklabels{l};
   else
      yl = yticklabels(l);
   end
   if yt < axy(4)-Dy*Ly && yt > axy(3),
      if isnumeric(yl),
         if isnan(yl),
            yl = '';
         elseif useformat == false,
            yl = num2str(yl);
         else
            yl = num2str(yl,ylformat);
         end
      end
      yll = length(yl);
      h = plot([Ypos Ypos+12/figwidth*Dx],[yt yt],plotoptions{:});
      set(h,'tag',tag)
      if abs(yt-Ypos) < 0.01*Dy, yt = Ypos+0.01*Dy; end
      h = text(Ypos-0.005*Dx,yt,yl,'HorizontalAlignment','right',...
               fontoptions{:});
      set(h,'tag',tag)
   end
end

axis(axy)
if ~holdon, hold off, end

% Copyright (c) Maarten Jansen
%
% This software is part of ThreshLab and is copyrighted material.
