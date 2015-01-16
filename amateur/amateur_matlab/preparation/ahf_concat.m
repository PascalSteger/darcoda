% ahf_concat.m: concat all previously prepared ahf_ variables

ahf_hostno(ahf_hostno>=0) = ahf_hostno(ahf_hostno>=0) + length(m_ahf_hostno);
m_ahf_hostno = cat(1,m_ahf_hostno,ahf_hostno);

m_ahf_offset = cat(1,m_ahf_offset,length(m_ahf_offset)*ones(length(ahf_hostno),1));

m_ahf_npart1 = cat(1,m_ahf_npart1,ahf_npart1);
m_ahf_n_gas = cat(1,m_ahf_n_gas, ahf_n_gas);
m_ahf_subno = cat(1,m_ahf_subno,ahf_subno+length(m_ahf_subno));
m_ahf_mvir = cat(1,m_ahf_mvir,ahf_mvir);
m_ahf_rtot = cat(1,m_ahf_rtot,ahf_rtot);

%m_cosjt = cat(1,m_cosjt,cosjt);
%m_ahf_cosjt = cat(1,m_ahf_cosjt,ahf_cosjt);
%m_ahf_costea = cat(1,m_ahf_costea,ahf_costea);
%m_ahf_costeb = cat(1,m_ahf_costeb,ahf_costeb);
%m_ahf_costec = cat(1,m_ahf_costec,ahf_costec);

m_ahf_cosrj = cat(1, m_ahf_cosrj, ahf_cosrj);
m_ahf_cpa = cat(1, m_ahf_cpa, ahf_cpa);
m_ahf_lambda = cat(1,m_ahf_lambda, ahf_lambda);
m_ahf_sigv = cat(1,m_ahf_sigv,ahf_sigv);

m_ahf_jx = cat(1, m_ahf_jx, ahf_jx);
m_ahf_jy = cat(1, m_ahf_jy, ahf_jy);
m_ahf_jz = cat(1, m_ahf_jz, ahf_jz);
