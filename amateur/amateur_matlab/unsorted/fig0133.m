nbins = 8;

idxhost = find(hostno==-2);
allmu = [];
for i=1:length(idxhost)
    idxsub = find( hostno==idxhost(i) );
    mu = mvir(idxsub)/mvir(idxhost(i));
    allmu = [allmu; mu];
end

xmin = log10(min(mvir)/min(mvir(idxhost)));
xspace = logspace(xmin,0,nbins+1);
[nbin,ibin] = histc(allmu,xspace);

mf  = [];
emf = [];
xm  = [];

for i=1:nbins
    ii = ibin==i;
    mf(i) = nbin(i)/(log(xspace(i+1))-log(xspace(i)));
    emf(i) = mf(i)/sqrt(nbin(i));
    xm(i) = median(allmu(ii));
end

errorbar(log10(xm),log10(mf),log10((mf-emf)./mf),log10((mf+emf)./mf));%,log10(emf));
brob = robustfit(log10(xm(3:end)),log10(mf(3:end)));
hold on;
xx = linspace(xmin,0,3);
plot(xx,brob(2)*xx+brob(1),'r-');
hold off;
%set(gca,'Xscale','log','Yscale','log');


