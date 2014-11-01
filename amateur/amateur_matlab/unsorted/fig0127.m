%amautinput
%amautprepare

% reproduction of fig. 4 Pereira et al. 2008
nbins=10;
rmax=3; %in r_{vir} of host halo
r=rtot; r2=r; mvir2=mvir;
%acp=abs(cosphirva0(r<rmax));
acp=abs(cosphirv3(r<rmax));

r2=r(r<rmax); mvir2=mvir(r<rmax);

%[nbin,ibin]=histc(rtot,linspace(0,max(rtot),nbins(j)+1));

%[nbin,ibin]=histc(r2,linspace(0,rmax,nbins+1));
[nbin,ibin]=histc(r2,prctile(r2,[linspace(0,100,nbins+1)]));

rmed1 = []; phimed1 = []; phierru1 = []; phierrl1 = [];
rmed2 = []; phimed2 = []; phierru2 = []; phierrl2 = [];
nnzj1 = []; nnzj2 = [];
for i=1:nbins
    %bnzj = npart1'>200*(j-1) & npart1'<200*j & ibin==i & ~isnan(acp);
    %bnzj = npart1'>2000*(j-1) & ibin==i & ~isnan(acp) & r2<1e5;
    
    mres=1e11;
    bnzj1 = mvir2'>10*min(mvir(npart1==min(npart1))) & mvir2'<mres & ibin==i & ~isnan(acp);
    nnzj1(i) = sum(bnzj1);
    
    bnzj2 = mvir2'>mres & ibin==i & ~isnan(acp);
    nnzj2(i) = sum(bnzj2);
    
    rmed1(i)=-1; phimed1(i)=-1;cos(0.0)*180/pi;
    rmed2(i)=-1; phimed2(i)=-1;
    rmed1(i) = median(r2(bnzj1));
    rmed2(i) = median(r2(bnzj2));
    phimed1(i) = mean(acp(bnzj1));
    phimed2(i) = mean(acp(bnzj2));
    phierru1(i) = prctile(acp(bnzj1),84)-phimed1(i);
    phierru2(i) = prctile(acp(bnzj2),84)-phimed2(i);
    phierrl1(i) = phimed1(i) - prctile(acp(bnzj1),16);
    phierrl2(i) = phimed2(i) - prctile(acp(bnzj2),16);
end
subplot(1,2,1)
errorbar(rmed1,phimed1,phierru1,phierrl1,'Marker','hexagram','LineStyle','-','Color',[1 0 0]);
hold on;
errorbar(rmed1,phimed1,(phierru1-phierrl1)/2./sqrt(nnzj1),'Marker','hexagram','LineStyle','-','Color',[1 0 1]);
scatter(r2(mvir2'>10*min(mvir(npart1==min(npart1))) & mvir2'<mres),acp(mvir2'>10*min(mvir(npart1==min(npart1))) & mvir2'<mres));
axis([0,rmax,0,1]);
 plot(0.5*ones(rmax,1),'--'); hold off;
xlabel('r / r_{vir}'); ylabel(sprintf('<|cos(r,L)|> @ m_{vir}< %d',mres));

subplot(1,2,2)
errorbar(rmed2,phimed2,phierru2,phierrl2,'Marker','hexagram','LineStyle','-','Color',[1 0 0]);
hold on;
errorbar(rmed2,phimed2,(phierru2-phierrl2)/2./sqrt(nnzj2),'Marker','hexagram','LineStyle','-','Color',[1 0 1]);
axis([0,rmax,0,1]);
scatter(r2(mvir2'>mres),acp(mvir2'>mres));
plot(0.5*ones(rmax,1),'--'); hold off;
xlabel('r / r_{vir}'); ylabel(sprintf('<|cos(r,L)|> @ m_{vir}> %d',mres));
