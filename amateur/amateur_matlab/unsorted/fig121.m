%start;
%pv_dm_gas = cosphiauto(1,:,10);pv_gas_stars = cosphiauto(1,:,9);pv_dm_stars = cosphiauto(1,:,16);

nbins=8;
rmax=3;
r=rtot; r2=r(r<rmax);

acp  = abs(cosphirea(1,:,4));
excm = exc(r<rmax);
mvir2= ahf_mvir(r<rmax); mres=1e11;


[nbin,ibin]=histc(r2,prctile(r2,[linspace(0,100,nbins+1)]));
rmed1 = []; phimed1 = []; phierru1 = []; phierrl1 = [];
rmed2 = []; phimed2 = []; phierru2 = []; phierrl2 = [];
nnzj1 = []; nnzj2 = [];


for i=1:nbins    
    bnzj1 = excm & mvir2<mres & ibin==i & ~isnan(acp(r<rmax));
    nnzj1(i) = sum(bnzj1);
    
    bnzj2 = excm & mvir2>mres & ibin==i & ~isnan(acp(r<rmax));
    nnzj2(i) = sum(bnzj2);
    
    rmed1(i)=-1; phimed1(i)=-1;
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
scatter(r2(excm & mvir2<mres),acp(excm & mvir2<mres));
axis([0,rmax,0,1]);
plot(0.5*ones(rmax,1),'--'); hold off;
xlabel('r/r_{vir}');
ylabel(sprintf('<|cos(r,L)|> @ m_{vir}< %d',mres));



subplot(1,2,2)

errorbar(rmed2,phimed2,phierru2,phierrl2,'Marker','hexagram','LineStyle','-','Color',[1 0 0]);
hold on; axis([0,rmax,0,1]); scatter(r2(mvir2>mres),acp(mvir2>mres)); plot(0.5*ones(rmax,1),'--'); hold off;
xlabel('r/r_{vir}'); ylabel(sprintf('<|cos(r,L)|> @ m_{vir}> %d',mres));
