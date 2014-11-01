% equivalent to fig. 3 of van den Bosch et al. 2002
% or better fig. 10 Croft et al. 2008
[n,x] = hist(acos((m_ahf_cpa(exc)))*180/pi,12);

bar(x,n/length((m_ahf_cpa(exc))),1);
xlabel('phi(L_{gas}, L_{DM})'); ylabel('fraction');

mean(acos((m_ahf_cpa(m_exc & ~isnan(m_ahf_cpa))))*180/pi)
median(acos((m_ahf_cpa(m_exc & ~isnan(m_ahf_cpa))))*180/pi)

x=[x', n'];
save '/data/achtland1/psteger/amd/halos/vis/report/fig110.dat' x '-ASCII';

%8