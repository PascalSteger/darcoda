% hp_lambdaprime.m
nbin = 20;
figprep('\lambda',...
        'log_{10}\lambda',...
        'log_{10}p(\lambda)');


exl = m_exc & m_hp_lambda(:,1,1)>0;
pl  = m_hp_lambda(exl,1,1);
hist(log10(pl),nbin);