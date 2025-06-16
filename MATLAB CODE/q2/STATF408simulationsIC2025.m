if exist('realrandom') ~= 1, realrandom = false; end
if exist('samplesize') ~= 1, samplesize = 30; end
if exist('fullmodelsize') ~= 1, fullmodelsize = 30; end
if exist('knots') ~= 1,
   knots = [-0.1000 0.1555 0.3143 0.5469 0.6903 0.8730 1.1];
end
truemodelsize = length(knots)+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('makeprint') ~= 1, makeprint = false; end
if exist('figsize') ~= 1, figsize = [1280 420]; end
if exist('figpaperpos') ~= 1, figpaperpos = [0 0 figsize/96]; end
forestgreen = [34 139 34]/255; fg = forestgreen;
dotoptg = {'o','markersize',2,'markeredgecolor',fg,'markerfacecolor',fg};
dotoptr = {'o','markersize',2,'markeredgecolor','r','markerfacecolor','r'};
linedotopt = {'-o','linewidth',1,'markersize',2,'color'};
linedotoptg = {linedotopt{:},fg,'markerfacecolor','w'};
linedotoptr = {linedotopt{:},'r','markerfacecolor','w'};
linedotoptb = {linedotopt{:},'b','markerfacecolor','w'};
lineoptb = {'linewidth',1,'color','b'};
lineoptr = {'linewidth',1,'color','r'};
lineoptg = {'linewidth',1,'color',fg};
fontsize = 10;
textopt = {'fontsize',fontsize,'fontweight','bold'};

if realrandom == false, rand('state',0); randn('state',0); end
n = samplesize;
m = fullmodelsize;
k = truemodelsize-1;  % k is # zeros of the polynomial: these k parameters fix
                      % the polynomial up to a constant 
makex = ~exist('x');
if exist('x') == 1, makex = length(x)~=n; end
if makex, x = sort(rand(n,1)); end
makeZ = ~exist('Z');
if exist('Z') == 1, makeZ = length(Z)~=n; end
if makeZ, Z = randn(n,1); end
% sometimes (in Octave), empirical variance of Z seems to be far from 1,
% leading to a Cp curve far from the PE curves. In that case, apply the
% following rescaling:
% Z = Z/sqrt(mean((Z-mean(Z)).^2));
if exist('stdev') ~= 1, stdev = 0.2; end

nn = 1001; xx = (0:nn-1)'/(nn-1);

mu = ones(size(x)); mmu = ones(size(xx));
for i = 1:k, mu = mu.*(x-knots(i)); mmu = mmu.*(xx-knots(i)); end
mmumax = max(abs(mmu));
mmu = mmu/mmumax;
mu = mu/mmumax;

Y = mu+stdev*Z;

figure(1)
plot(x,Y,dotoptg{:})
hold on
plot(xx,mmu,lineoptr{:})
plot(xx,zeros(size(xx)),'k')
hold off
axis off
figpos = get(gcf,'position'); figpos(3:4) = figsize;
set(gcf,'position',figpos,'paperposition',figpaperpos)
drawaxes('xtick',[0 1],'ytick',[-1 1],textopt{:})

if makeprint == true,
   print -dpdf nestedmodel.pdf
end


% PE and Cp are defined without studentization or standardization
% i.e.: PE = E(norm(mu-muhat)^2)/n
%  and: Cp = SSE/n+2p/n*stdev^2-stdev^2
X = ones(n,m); models = cell(1,m); models{1} = 1;
for p=2:m, X(1:n,p) = (x-0.5).^(p-1); models{p} = (1:p); end
pp = (1:m);

PE = PElinearmodel(mu,X,models,stdev);
[Cp AIC BIC betahats] = CpAICBIClinearmodel(Y,X,models,stdev);
Cp = Cp*stdev^2/n;
[dif2 KL] = squareddifKLlinearmodel(Y,X,models,mu,stdev);
% squared difference is norm(muhat-mu)
% PE = Expected squared difference

figure(2)
% plots squared difference, PE and Cp as function of p
plot(pp,PE,linedotoptg{:})
hold on
plot(pp,dif2,linedotoptb{:})
plot(pp,Cp,linedotoptr{:})
hold off
aa = axis; % aa(3) = 0;
axis(aa)

xticks = get(gca,'XTick');
yticks = get(gca,'YTick');

axis off
aa = axis; aa(2) = 1.1*aa(2); aa(4) = 1.2*aa(4)-0.2*aa(3); axis(aa)
drawgrids(xticks,yticks,'k-')
drawaxes('linewidth',2,'xticks',xticks,textopt{:})
Dx = aa(2)-aa(1);
Dy = aa(4)-aa(3);
text(aa(2),-fontsize/150*Dy,'Model size (p)',...
          'HorizontalAlignment','right',textopt{:})
h = legend('PE','squared diff','C_p');
set(h,textopt{:})

figpos = get(gcf,'position'); figpos(3:4) = figsize;
set(gcf,'position',figpos,'paperposition',figpaperpos)
drawaxes('linewidth',2,'xticks',xticks,textopt{:})

if makeprint == true,
   print -dpdf Cpnestedmodel.pdf
end

% Now extend the plot of beyond [0,1]
figure(3)
xx = xx*(knots(end)-knots(1)) + knots(1);
[minPE popt] = min(PE);
% Xopt = X(1:n,1:popt); betahatopt = (Xopt'*Xopt)\(Xopt'*Y);
% muhatopt = Xopt*betahatopt;
betahatopt = betahats{popt};
XXopt = ones(nn,1);
for p=1:popt-1, XXopt = [XXopt (xx-0.5).^p]; end
mmuhatopt = XXopt*betahatopt;
mmu = ones(size(xx))/mmumax;
for i = 1:k, mmu = mmu.*(xx-knots(i)); end

plot(x,Y,dotoptg{:})
hold on
plot(xx,mmu,lineoptr{:})
plot(xx,mmuhatopt,lineoptb{:})
hold off
axis([knots(1) knots(end) -1.1 1.1])
axis off
figpos = get(gcf,'position'); figpos(3:4) = figsize;
set(gcf,'position',figpos,'paperposition',figpaperpos)
drawaxes('xtick',[0 1],'ytick',[-1 1],textopt{:})

if makeprint == true,
   print -dpdf minPEestnestedmodel.pdf
   close all
end

figure(4)
ElogLDGP = -1/2-log(sqrt(2*pi)*stdev);
plot(pp,KL-ElogLDGP,linedotoptb{:})
hold on
plot(pp,-AIC/2,linedotoptr{:})
aa = axis;
plot(pp,-BIC/2,linedotoptg{:})
hold off
axis(aa)

xticks = get(gca,'XTick');
yticks = get(gca,'YTick');

axis off
aa = axis; aa(2) = 1.1*aa(2); aa(4) = 1.2*aa(4); axis(aa)
drawgrids(xticks,yticks,'k-')
drawaxes('linewidth',2,'xticks',xticks,textopt{:})
Dx = aa(2)-aa(1);
Dy = aa(4)-aa(3);
text(aa(2),-fontsize/150*Dy,'Model size (p)',...
          'HorizontalAlignment','right',textopt{:})
h = legend('KL-ElogL(DGP)','-AIC','-BIC');
set(h,textopt{:})

figpos = get(gcf,'position'); figpos(3:4) = figsize;
set(gcf,'position',figpos,'paperposition',figpaperpos)
drawaxes('linewidth',2,'xticks',xticks,textopt{:})

