
n = 10000;
X = randexpinvsqrt(n,1);
h = 0.1; kerneltype = 'cos';
u = (0.5+(0:n-1))/n;
x = 2*u-1;
fXhat = kernelestimation(x,X,h,kerneltype);
figure(2)
f = exp(-1./(1-x.^2));
dx = initmom(x,0);
K = 1/(f*dx);
fX = K*f;
plot(x,fXhat,'b','linewidth',3)
hold on
plot(x,fX,'r','linewidth',3)
hold off
axy = axis; axy(2) = 1.1; axy(4) = axy(4)*1.1; axis(axy);
fticks = bestticks(0,axy(4)*1.1,10);
axis off
drawaxes
drawgrids(xticks,fticks)

varX = (fX.*x.^2)*dx
varXhat = mean(X.^2)

% The variance estimator may not be that accurate, but the procedure allows us
% to generate pseudo-random numbers.
