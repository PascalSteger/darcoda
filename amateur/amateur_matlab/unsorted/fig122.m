%start;
%pv_dm_gas = cosphiauto(1,:,10);pv_gas_stars = cosphiauto(1,:,9);pv_dm_stars = cosphiauto(1,:,16);

nbins=12;
r=m_rtot; rmax=3;

excr = r<rmax;
excs = m_ahf_hostno<=0;
%excs = m_ahf_hostno<Inf;

excp = excr & excs;
r2=r(excp);

acp  = abs(m_cosphirea(1,:,1));
%acp  = abs(m_cosphirec(1,:,4));
excm = m_exc(excp);
mvir2= m_hp_mvir(excp); mres=10;


[nbin,ibin]=histc(r2,prctile(r2,[linspace(0,100,nbins+1)]));
rmed1 = []; phimed1 = []; phierru1 = []; phierrl1 = [];
rmed2 = []; phimed2 = []; phierru2 = []; phierrl2 = [];
nnzj1 = []; nnzj2 = [];

acp2 = acp(excp);

for i=1:nbins    
    bnzj1 = excm & mvir2<mres & ibin==i & ~isnan(acp2);
    nnzj1(i) = sum(bnzj1);
    
    bnzj2 = excm & mvir2>mres & ibin==i & ~isnan(acp2);
    nnzj2(i) = sum(bnzj2);
    
    rmed1(i)=-1; phimed1(i)=-1;
     rmed2(i)=-1; phimed2(i)=-1;
     
    rmed1(i) = median(r2(bnzj1));
     rmed2(i) = median(r2(bnzj2));
     
    phimed1(i) = median(acp2(bnzj1));
     phimed2(i) = median(acp2(bnzj2));
     
    phierru1(i) = prctile(acp2(bnzj1),84)-phimed1(i);
     phierru2(i) = prctile(acp2(bnzj2),84)-phimed2(i);
    phierrl1(i) = phimed1(i) - prctile(acp2(bnzj1),16);
     phierrl2(i) = phimed2(i) - prctile(acp2(bnzj2),16);
end



subplot(1,2,1)
scatter(r2(excm & mvir2<mres),acp2(excm & mvir2<mres),'.'); hold on; plot(0.5*ones(rmax,1),'--');
errorbar(rmed1,phimed1,phierru1,phierrl1,'Marker','hexagram','LineStyle','-','Color',[1 0 0]);hold off;
axis([0,rmax,0,1]);
xlabel('r/r_{vir}');
ylabel(sprintf('<|cos(r,L)|> @ m_{vir}< %d',mres));


subplot(1,2,2)
scatter(r2(mvir2>mres),acp2(mvir2>mres),'.'); hold on; plot(0.5*ones(rmax,1),'--');
errorbar(rmed2,phimed2,phierru2,phierrl2,'Marker','hexagram','LineStyle','-','Color',[1 0 0]);hold off;
axis([0,rmax,0,1]);  
xlabel('r/r_{vir}'); ylabel(sprintf('<|cos(r,L)|> @ m_{vir}> %d',mres));
