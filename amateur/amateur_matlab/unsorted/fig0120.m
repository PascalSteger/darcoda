%amautprepare

% reproduction of fig. 4 Pereira et al. 2008
acp=abs(cosphirv0);

nbins = 10;
[nbin,ibin]=histc(rtot,linspace(0,max(rtot),nbins+1));
rmed = []; phimed = []; phierru = []; phierrl = [];
for i=1:nbins
    j = ibin==i & ~isnan(acp);
    rmed(i) = median(rtot(j));
    phimed(i) = median(acp(j));
    phierru(i) = prctile(acp(j),84)-phimed(i);
    phierrl(i) = phimed(i) - prctile(acp(j),16);
end
errorbar(rmed,phimed,phierru, phierrl,'s')