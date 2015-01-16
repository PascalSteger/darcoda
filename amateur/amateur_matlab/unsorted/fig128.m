%start;
%pv_dm_gas = cosphiauto(1,:,10);pv_gas_stars = cosphiauto(1,:,9);pv_dm_stars = cosphiauto(1,:,16);

nbins=8;
excr = m_rtot<30;
%excs = m_ahf_hostno<=0;
excs = m_ahf_hostno<Inf;

excp = excr & excs;
xax  = m_Tew(1,excp)-1/3*sum(m_Tew(:,excp),1);
%xax  = m_Tew(3,excp);
xmax = max(xax);
xmin = min(xax)/10;

acp  = abs(m_cosphilt(1,:,1)); %abc,:,1-6

%acp  = acp(exc);
excm = m_exc(excp);
mvir2= m_ahf_mvir(excp); mres=1e11;


[nbin,ibin]=histc(xax,prctile(xax,[linspace(0,100,nbins+1)]));
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
     
    rmed1(i) = median(xax(bnzj1));
     rmed2(i) = median(xax(bnzj2));
     
    phimed1(i) = median(acp2(bnzj1));
     phimed2(i) = median(acp2(bnzj2));
     
    phierru1(i) = prctile(acp2(bnzj1),84)-phimed1(i);
     phierru2(i) = prctile(acp2(bnzj2),84)-phimed2(i);
    phierrl1(i) = phimed1(i) - prctile(acp2(bnzj1),16);
     phierrl2(i) = phimed2(i) - prctile(acp2(bnzj2),16);
end



subplot(1,2,1)
scatter(xax(excm & mvir2<mres),acp2(excm & mvir2<mres),1,'.'); hold on; plot(0.5*ones(rmax,1),'--');
%errorbar(rmed1,phimed1,phierru1,phierrl1,'Marker','hexagram','LineStyle','-','Color',[1 0 0]);
errorbar(rmed1,phimed1,phierru1./sqrt(nnzj1),phierrl1./sqrt(nnzj1),'Marker','hexagram','LineStyle','-','Color',[1 0 0]);hold off;
axis([xmin,xmax,0,1]);
xlabel('T1');
ylabel(sprintf('<|cos(r,L)|> @ m_{vir}< %d',mres));


subplot(1,2,2)
scatter(xax(mvir2>mres),acp2(mvir2>mres),1,'.'); hold on; plot(0.5*ones(rmax,1),'--');
%errorbar(rmed2,phimed2,phierru2,phierrl2,'Marker','hexagram','LineStyle','-','Color',[1 0 0]);
errorbar(rmed2,phimed2,phierru2./sqrt(nnzj2),phierrl2./sqrt(nnzj2),'Marker','hexagram','LineStyle','-','Color',[1 0 0]);hold off;
axis([xmin,xmax,0,1]);  
xlabel('T1'); ylabel(sprintf('<|cos(r,L)|> @ m_{vir}> %d',mres));
