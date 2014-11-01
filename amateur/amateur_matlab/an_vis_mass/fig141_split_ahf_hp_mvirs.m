clf
set(gca,'FontSize',15);

subplot(1,2,1);
loglog(m_hp_mvir,m_ahf_mvir,'.')

subplot(1,2,2);
exmvir = m_ahf_mvir./m_hp_mvir>1.1e10;
loglog(m_hp_mvir(exmvir),m_ahf_mvir(exmvir),'.');
