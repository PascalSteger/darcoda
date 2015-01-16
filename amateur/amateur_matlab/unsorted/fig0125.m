%amautinput
%amautprepare

% reproduction of fig. 4 Pereira et al. 2008
nbins = 15;
rmax = 3; %in r_{vir} of host halo
r=rtotold; r2 = r; mvir2 = ahf_mvir;
cprl = cosphirl(1,exc' & r<rmax,1);
cprl = cosphirv(1,exc' & r<rmax,1);
acp = abs(cprl);
r2 = r(exc' & r<rmax); mvir2 = ahf_mvir(exc' & r<rmax);

[nbin,ibin] = histc(r2,linspace(0,rmax,nbins+1));
rmed1 = []; phimed1 = []; phierru1 = []; phierrl1 = [];
rmed2 = []; phimed2 = []; phierru2 = []; phierrl2 = [];
for i=1:nbins
    mdiff=1e12;
    bnzj1 = mvir2'<mdiff & ibin==i & ~isnan(acp);
    sm1 = sum(bnzj1);
    bnzj2 = mvir2'>mdiff & ibin==i & ~isnan(acp);
    sm2 = sum(bnzj2);
    
    rmed1(i)=-1; phimed1(i)=-1;
    rmed2(i)=-1; phimed2(i)=-1;
    rmed1(i) = median(r2(bnzj1));
    rmed2(i) = median(r2(bnzj2));
    phimed1(i) = median(acp(bnzj1));
    phimed2(i) = median(acp(bnzj2));
    phierru1(i) = (prctile(acp(bnzj1),84)-phimed1(i))/sqrt(sm1);
    phierru2(i) = (prctile(acp(bnzj2),84)-phimed2(i))/sqrt(sm2);
    phierrl1(i) = (phimed1(i) - prctile(acp(bnzj1),16))/sqrt(sm1);
    phierrl2(i) = (phimed2(i) - prctile(acp(bnzj2),16))/sqrt(sm2);
end
subplot(1,2,1)
errorbar(rmed1,phimed1,phierru1,phierrl1,'Marker','hexagram','LineStyle','-','Color',[1 0 0]);
axis([0,rmax,0,1]);
hold on; plot(0.5*ones(rmax,1),'--'); hold off;
xlabel('r / r_{vir}'); ylabel(sprintf('<|cos(r,L)|> @ m_{vir}< %d',mdiff));

subplot(1,2,2)
errorbar(rmed2,phimed2,phierru2,phierrl2,'Marker','hexagram','LineStyle','-','Color',[1 0 0]);
axis([0,rmax,0,1]);
hold on; plot(0.5*ones(rmax,1),'--'); hold off;
xlabel('r / r_{vir}'); ylabel(sprintf('<|cos(r,L)|> @ m_{vir}> %d',mdiff));
