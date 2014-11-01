% hp_concat.m: concat all previously prepared hp_ variables
% requires that the corresponding m_ arrays 
% were already declared in prepare_arrays.m

m_hp_mtot = cat(1,m_hp_mtot,hp_mtot);
m_hp_mvir = cat(1,m_hp_mvir,hp_mvir);
m_hp_npart= cat(1,m_hp_npart,hp_npart);
m_hp_dhost= cat(1,m_hp_dhost,hp_dhost);
m_hp_dhosttot = cat(1,m_hp_dhosttot, hp_dhosttot);
m_hp_dhostvir = cat(1,m_hp_dhostvir, hp_dhostvir);
m_hp_dhostvirtot = cat(1,m_hp_dhostvirtot, hp_dhostvirtot);
m_hp_dbig = cat(1,m_hp_dbig, hp_dbig);
m_hp_dbigtot = cat(1,m_hp_dbigtot, hp_dbigtot);

m_hp_j = cat(1,m_hp_j, hp_j);
m_hp_jtot = cat(1,m_hp_jtot, hp_jtot);
m_hp_lambda = cat(1,m_hp_lambda, hp_lambda);
m_hp_sigma  = cat(1,m_hp_sigma,  hp_sigma );

m_hp_eatot  = cat(1,m_hp_eatot,  hp_eatot);
m_hp_ebtot  = cat(1,m_hp_ebtot,  hp_ebtot);
m_hp_ectot  = cat(1,m_hp_ectot,  hp_ectot);