start;

[n,x] = hist(abs(ahf_cpa(exc)),10);

bar(x,n/sum(~isnan(ahf_cpa(exc))),1);

xlabel('|cos(L_{gas}, L_{DM})|'); ylabel('fraction');