% equivalent to fig. 3 of van den Bosch et al. 2002
% or better fig. 10 Croft et al. 2008

%start;

hp_cpa = cosphiauto(1,:,2);

[n,x] = hist(acos((hp_cpa(exc)))*180/pi,20);

bar(x,n/length((hp_cpa(exc))),1);
bar(x,n,1);

xlabel('phi(L_{gas}, L_{DM})'); ylabel('fraction');

mean(acos((hp_cpa(exc & ~isnan(hp_cpa))))*180/pi)
median(acos((hp_cpa(exc & ~isnan(hp_cpa))))*180/pi)