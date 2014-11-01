% tf_concat.m: concat all previously prepared ?[tf_,*]variables
% requires that the corresponding m_ arrays were already declared

m_hp_cosauto = cat(1,m_hp_cosauto, hp_cosauto);

m_hp_cosrj = cat(1,m_hp_cosrj, hp_cosrj);
m_hp_cosjea = cat(1,m_hp_cosjea, hp_cosjea);
m_hp_cosrea = cat(1,m_hp_cosrea, hp_cosrea);
m_hp_cosreb = cat(1,m_hp_cosreb, hp_cosreb);
m_hp_cosrec = cat(1,m_hp_cosrec, hp_cosrec);

m_hp_dhosttot = cat(1,m_hp_dhosttot,hp_dhosttot);
m_hp_cosjt = cat(1,m_hp_cosjt,hp_cosjt);
m_hp_cosrt = cat(1,m_hp_cosrt,hp_cosrt);

m_Tew = cat(1,m_Tew, Tew);
