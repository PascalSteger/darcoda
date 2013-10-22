#!/usr/bin/python

x = [trnd(20,1,50) trnd(4,1,100)+3]
hist(x,-2.25:.5:7.25)

pdf_normmixture = @(x,p,mu1,mu2,sigma1,sigma2) p*normpdf(x,mu1,sigma1) + (1-p)*normpdf(x,mu2,sigma2)


pStart = .5
muStart = quantile(x,[.25 .75])
sigmaStart = sqrt(var(x) - .25*diff(muStart).^2)
start = [pStart muStart sigmaStart sigmaStart]

lb = [0 -Inf -Inf 0 0]
ub = [1 Inf Inf Inf Inf]

paramEsts = mle(x, 'pdf',pdf_normmixture, 'start',start, 'lower',lb, 'upper',ub)

options = statset('MaxIter',300, 'MaxFunEvals',600)
paramEsts = mle(x, 'pdf',pdf_normmixture, 'start',start,'lower',lb, 'upper',ub, 'options',options)

bins = -2.5:.5:7.5
h = bar(bins,histc(x,bins)/(length(x)*.5),'histc')
set(h,'FaceColor',[.9 .9 .9])
xgrid = linspace(1.1*min(x),1.1*max(x),200)
pdfgrid = pdf_normmixture(xgrid,paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5))
hold on
plot(xgrid,pdfgrid,'-')
hold off
xlabel('x')
ylabel('Probability Density')
