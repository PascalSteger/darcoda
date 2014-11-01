% fig132_cosauto_asfct_r.m: cosauto as fct of radius
clf
set(gca,'FontSize',15);

m_cosauto_dm_gas = m_cosauto(1,:,7);
m_cosauto_gas_stars = m_cosauto(1,:,8);
m_cosauto_dm_stars = m_cosauto(1,:,12);

% KNOB
nbins=6;
r=m_rtot;
% KNOB
rmax=4;

excr = r<rmax;
% KNOB
%excs = m_ahf_hostno==-2;
excs = m_ahf_hostno<Inf;

excp = excr & excs;
r2=r(excp);

% KNOB
acp = abs(m_cosauto_dm_gas);
%acp = abs(m_cosauto_gas_stars);
acp = abs(m_cosauto_dm_stars);

% change as well.
excm = m_exc_4(excp);
mvir2= m_hp_mvir(excp);

% KNOB
mres=10; % x 10^{10}M_o


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
errorbar(rmed1,phimed1,phierrl1,phierru1,'Marker','hexagram','LineStyle','-','Color',[1 0 0]);hold off;
axis([0,rmax,0,1]);
xlabel('r/r_{vir}');
ylabel(sprintf('<|cos(j_k,j_l)|> @ m_{vir}< %d x 10^{10}M_o',mres));
grid on;

subplot(1,2,2)
scatter(r2(mvir2>mres),acp2(mvir2>mres),'.'); hold on; plot(0.5*ones(rmax,1),'--');
errorbar(rmed2,phimed2,phierrl2,phierru2,'Marker','hexagram','LineStyle','-','Color',[1 0 0]);hold off;
axis([0,rmax,0,1]);  
xlabel('r/r_{vir}');
ylabel(sprintf('<|cos(j_k,j_l)|> @ m_{vir}> %d x 10^{10}M_o',mres));
grid on;

%hist(acp,20)