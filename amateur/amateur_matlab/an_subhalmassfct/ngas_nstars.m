figprep('overview: number of star and gas particles. red: excluded',...
        'log_{10}n_{gas}','log_{10}n_{star}');
ng = m_hp_npart(:,1,3);
ns = m_hp_npart(:,1,4);
scatter(log10(ng),log10(ns),'.r')

exl = m_exc_3 & m_exc_4;
scatter(log10(ng(exl)),log10(ns(exl)),'.b')