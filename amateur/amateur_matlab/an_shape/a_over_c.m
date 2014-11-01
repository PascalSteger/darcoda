%a_over_c.m: plot ahf and hp histlogs for a/c
nbins = 20;
figprep('a/c','log_{10}a/c','log_{10}p(a/c)');

hold on;
histlog(ahf_a./ahf_c,nbins,'b');
histlog(m_hp_eatot(m_exc,1,1)./m_hp_ectot(m_exc,1,1),nbins,'r')
%[n,xout] = hist(log10(hp_eatot(:,1,1)./hp_ectot(:,1,1)),nbins);
%plot(xout,log10(n),'r');
%errorbar(xout,log10(n),log10(n)./sqrt(n),'ro');

legend('ahf','','hp','');