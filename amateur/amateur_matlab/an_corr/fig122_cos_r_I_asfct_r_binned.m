% 122_cor_r_I_asfct_r_binned.m: radial dependence of correlation between angles between r and I, binned in two mass bins
%pv_dm_gas = cosauto(1,:,10);pv_gas_stars = cosauto(1,:,9);pv_dm_stars = cosauto(1,:,16);
clf

% KNOB
k=2; %halotype
acp = abs(m_cosrea(1,:,k)); % 1,:,halotype(2: gas, 3: dm, 4: stars) % signal! use nbins=9 for 3, nbins=6 else
%acp = abs(m_cosrt(1,:,1)); % 1,:,Teabc; correlation between maximum tidal
%field eigenvector and direction vector towards main halo
%acp = abs(m_cosreb(1,:,1)); % signal for heavy halos
%acp = abs(m_cosrec(1,:,1)); % signal!
%acp = abs(m_cosjt(1,:,k)); %123,:,halotype
%acp = abs(m_cosrj(1,:,k));  %1,:,halotype

% KNOB
nbins=8;
r=m_rtot;
% KNOB
rmax=3; %*rvir

excr = r<rmax;
% KNOB
excs = m_ahf_hostno > -2; %hosts
%excs = m_ahf_hostno ==-1; %single/field
%excs = m_ahf_hostno >= 0; %subs
% take all halos and subhalos from AHFstep:
%excs = ones(size(m_ahf_hostno));

if k==1
    exct = m_exc_1;
else if k==2
        exct = m_exc_2;
    else if k == 3
            exct = m_exc_3;
        else if k == 4
                exct = m_exc_4;
            end
        end
    end
end

excp = excr & excs & exct;
%excp = excr & excs;

r2=r(excp);
excm = m_exc_1(excp);
mvir2= m_hp_mvir(excp);

% KNOB
mres=10000000; % x 10^{10}M_o

[nbin,ibin]=histc(r2,prctile(r2,[linspace(0,100,nbins+1)]));
rmed1 = []; phimed1 = []; phierru1 = []; phierrl1 = [];
rmed2 = []; phimed2 = []; phierru2 = []; phierrl2 = [];
nnzj1 = []; nnzj2 = [];

acp2 = acp(excp);

for i=1:nbins    
    bnzj1 = excm & mvir2 < mres & ~isnan(acp2) & ibin==i;
    nnzj1(i) = sum(bnzj1);
    
    bnzj2 = excm & mvir2 > mres & ~isnan(acp2) & ibin==i;
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

%subplot(1,2,1);set(gca,'FontSize',15);
scatter(r2(excm & mvir2<mres),acp2(excm & mvir2<mres),'.'); hold on;
errorbar(rmed1,phimed1,phierrl1,phierru1,'Marker','o','LineStyle','-','Color',[1 0 0]);hold off;
axis([0,rmax,0,1]);
xlabel('r/r_{vir}');
ylabel(sprintf('<|cos(r,e_a)|> @ m_{vir}< %d x 10^{10}M_o',mres));
grid on;

%subplot(1,2,2);set(gca,'FontSize',15);
%scatter(r2(mvir2>mres),acp2(mvir2>mres),'.'); hold on;
%errorbar(rmed2,phimed2,phierrl2,phierru2,'Marker','o','LineStyle','-','Color',[1 0 0]);hold off;
%axis([0,rmax,0,1]);  
%xlabel('r/r_{vir}');
%ylabel(sprintf('<|cos(r,e_a)|> @ m_{vir}> %d x 10^{10}M_o',mres));
%grid on;

%% searching for other common properties of halos with cos angle in upper/lower third

%excpo1 = acp2(logical(excm))>0.66;
%excpo2 = acp2(logical(excm))<0.33;
%subplot(1,2,1); hist(r2(excpo1),10);
%subplot(1,2,2); hist(r2(excpo2),10);
