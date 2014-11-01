% 110_crosscorr_jgasdm.m: distribution of cos(angle) between angular momenta of gas and dm components in a halo
% equivalent to fig. 3 of van den Bosch et al. 2002
% or more recently fig. 10 Croft et al. 2008
clf;
set(gca,'FontSize',15);

cpa = m_ahf_cpa;
%cpa = m_hp_cpa;
cut1 = 1; cut2 = 3;

% xor split for distance bins
%exl1 = m_exc_1 & m_rtot < cut1;
%exl2 = m_exc_1 & m_rtot > cut1 & m_rtot < cut2;
%exl3 = m_exc_1 & m_rtot > cut2;
% xor according to halo type
mth = [];
for k=1:length(m_ahf_hostno);
    if m_ahf_hostno(k) < 0
        mth(k) = m_ahf_hostno(k);
        continue;
    end
    if m_ahf_hostno(m_ahf_hostno(k)+1) < 0
        mth(k) = m_ahf_hostno(k);
    else
        mth(k) = -3;
    end
end
exl1 = m_exc_1 & m_ahf_n_gas > 0 & mth == -2;
exl2 = m_exc_1 & m_ahf_n_gas > 0 & mth  > -2;
% xor depending on mass
%msplit = prctile(m_ahf_mvir,[1/3 2/3]*100);
%m1 = msplit(1); m2 = msplit(2);
%exl1 = m_exc_1 & m_ahf_n_gas > 0 & m_ahf_mvir <= m1;
%exl2 = m_exc_1 & m_ahf_n_gas > 0 & m_ahf_mvir > m1 & m_ahf_mvir <= m2;
%exl3 = m_exc_1 & m_ahf_n_gas > 0 & m_ahf_mvir > m2;

cpa1 = cpa(exl1);
cpa2 = cpa(exl2);

nbin1 = 6;
nbin2 = 6;

[n1,x1] = hist(acos(cpa1)*180/pi,nbin1);
[n2,x2] = hist(acos(cpa2)*180/pi,nbin2);

subplot(1,2,1);set(gca,'FontSize',15);

bar(x1,n1/sum(n1),1);
%xlabel(sprintf('\phi (J_{gas}, J_{DM}) for r_{vir} < %d',cut1));
xlabel('\phi (J_{gas}, J_{DM}) for host halos');
ylabel('frequency');
%legend(sprintf('median: %d',
median(acos(cpa1)*180/pi)
%));

grid on;
axis([0, 180, 0, 0.3]);

subplot(1,2,2);set(gca,'FontSize',15);
bar(x2,n2/sum(n2),1);
%xlabel(sprintf('\phi (J_{gas}, J_{DM}) for %d < r_{vir} < %d',cut1,cut2));
xlabel('\phi (J_{gas}, J_{DM}) for single halos');
ylabel('frequency');
grid on;
%legend(sprintf('median: %d',
median(acos(cpa2)*180/pi)
%));
axis([0, 180, 0, 0.3]);

%mean(acos((m_ahf_cpa(m_exc_1 & ~isnan(m_ahf_cpa))))*180/pi)
%median(acos((m_ahf_cpa(m_exc_1 & ~isnan(m_ahf_cpa))))*180/pi)


% output for SuperMongo
%x=[x', n'];
%save '/data/achtland1/psteger/amd/halos/vis/report/fig110.dat' x '-ASCII'
