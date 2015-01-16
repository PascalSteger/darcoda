hh = double(hp_ptypn(:,1,3)) + double(hp_ptypn(:,1,4)) + double(hp_ptypn(:,1,6));
hh = hh./double(hp_ptypn(:,1,2));
hist(hh(~isnan(hh) & hh>0 & hh<1e-1),50);
xlabel('#lowres/#DM');
ylabel('a.u.');