studnumber = 12345;
if exist('realrandom') ~= 1, realrandom = false; end
if exist('samplesize') ~= 1, samplesize = 400; end
if exist('fullmodelsize') ~= 1, fullmodelsize = 2000; end
if exist('bandwidth') ~= 1, bandwidth = 0.1; end
if exist('degreeofsparsity') ~= 1, degreeofsparsity = 0.05; end
if exist('setupdesign') ~= 1, setupdesign = true; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Technical stuff
if exist('figsize') ~= 1, figsize = [1280 420]; end
if exist('figpaperpos') ~= 1, figpaperpos = [0 0 figsize/96]; end
forestgreen = [34 139 34]/255; fg = forestgreen;
dotoptg = {'o','markersize',2,'markeredgecolor',fg,'markerfacecolor',fg};
dotoptr = {'o','markersize',2,'markeredgecolor','r','markerfacecolor','r'};
linedotopt = {'-o','linewidth',1,'markersize',2,'color'};
linedotoptg = {linedotopt{:},fg,'markerfacecolor','w'};
linedotoptr = {linedotopt{:},'r','markerfacecolor','w'};
linedotoptb = {linedotopt{:},'b','markerfacecolor','w'};
linedotoptm = {linedotopt{:},'m','markerfacecolor','w'};
linedotoptk = {linedotopt{:},'k','markerfacecolor','w'};
lineoptb = {'linewidth',1,'color','b'};
lineoptr = {'linewidth',1,'color','r'};
lineoptg = {'linewidth',1,'color',fg};
fontsize = 10;
textopt = {'fontsize',fontsize,'fontweight','bold'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if realrandom==false, rand('state',studnumber); randn('state',studnumber); end
n = samplesize;
m = fullmodelsize;
p = degreeofsparsity;
stdev = 1;
if setupdesign | exist('X') ~= 1,
   x = sort(rand(n,1));
   u = sort(rand(1,m));
   h = bandwidth;
   xu = x*ones(1,m)-ones(n,1)*u;
   X = evalkernel(xu,h,'cos');
   colnormX = sqrt(sum(X.^2,1)); X = X*diag(1./colnormX);
end

image(X*255)
colormap gray
beta = zeros(m,1);
I = find(rand(size(beta))<p);
kappa = length(I);
beta(I) = log(rand(size(I)))*5.*sign(rand(size(I))-0.5);
plot(beta) % pause here to see beta plot
mu = X*beta;
Z = randn(size(mu));
Y = mu+stdev*Z;

[betahat muhat S Cp1 lambda Sopt] = LARS(Y,X,stdev);

nS = length(S)+1;
models = cell(1,nS);
kappa = zeros(1,nS);
modelk = []; models{1} = modelk;
for k = 1:nS-1,
   Sk = S(k);
   if Sk<0, modelk = setdiff(modelk,-Sk); else, modelk = [modelk Sk]; end
   models{k+1} = modelk;
   kappa(k+1) = length(modelk);
end
Cp0 = CpAICBIClinearmodel(Y,X,models,stdev)*stdev^2/n;
PE0 = PElinearmodel(mu,X,models,stdev);
samplePE0 = squareddifKLlinearmodel(Y,X,models,mu,stdev);
samplePE1 = squareddifLasso(Y,X,models,lambda,mu,stdev);
figure(1)
plot(kappa,samplePE0,linedotoptr{:})
hold on
plot(kappa,samplePE1,linedotoptm{:})
plot(kappa,Cp1,linedotoptb{:})
plot(kappa,Cp0,linedotoptk{:})
hold off
aa = axis; % aa(3) = 0;
axis(aa)
xticks = get(gca,'XTick');
yticks = get(gca,'YTick');

axis off
aa = axis; aa(2) = 1.1*aa(2); aa(4) = 1.2*aa(4)-0.2*aa(3); axis(aa)
drawgrids(xticks,yticks,'k-')
drawaxes('linewidth',1,'xticks',xticks,textopt{:})
Dy = aa(4)-aa(3);
text(aa(2),-fontsize/150*Dy,'Model size (\kappa)',...
          'HorizontalAlignment','right',textopt{:})
h = legend('PE orth.proj','PE shrinkage','C_p shrinkage','C_p orth.proj');
set(h,textopt{:})

figpos = get(gcf,'position'); figpos(3:4) = figsize;
set(gcf,'position',figpos,'paperposition',figpaperpos)
drawaxes('linewidth',1,'xticks',xticks,textopt{:})
