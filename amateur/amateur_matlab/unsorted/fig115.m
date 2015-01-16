%plot angle between different L's from hp

% make two or three mass bins
%exd = 
nbins = 50;

% fast access to different angle correlations
pv_dm_gas = (m_cosphiauto(1,:,8));
pv_gas_stars = abs(m_cosphiauto(1,:,7));
pv_dm_stars = abs(m_cosphiauto(1,:,12));

xax  = m_Tew(3,:)-1/3*sum(m_Tew(:,:),1);

splitval = sum(m_Tew(:,:),1);

valid = ~isnan(pv_dm_stars)& ~isnan(splitval) & m_ahf_hostno<0;
dbin1 = splitval<prctile(splitval(valid),33) & valid;
dbin2 = splitval>prctile(splitval(valid),33) &splitval<prctile(splitval(valid),67) & valid;
dbin3 = splitval>prctile(splitval(valid),67) & valid;

figure
hold on;
[n,x] = hist(pv_dm_gas(dbin1),nbins);stairs(x,n./length(~isnan(pv_dm_gas(dbin1))),'b-');
[n,x] = hist(pv_dm_gas(dbin2),nbins);stairs(x,n./length(~isnan(pv_dm_gas(dbin2))),'r-');

figure
hold on;
h=cdfplot(pv_dm_gas(dbin1)); set(h,'color','b');
h=cdfplot(pv_dm_gas(dbin2)); set(h,'color','g');
h=cdfplot(pv_dm_gas(dbin3)); set(h,'color','r');

figure;
[f1,x1] = ksdensity(pv_dm_gas(dbin1),'support',[-1 1]);
[f2,x2] = ksdensity(pv_dm_gas(dbin2),'support',[-1 1]);
[f3,x3] = ksdensity(pv_dm_gas(dbin3),'support',[-1 1]);
plot(x1,f1,x2,f2,x3,f3);
legend('1','2','3');

%figure
%[n,x] = hist(pv_gas_stars,nbins);bar(x,n./length(~isnan(pv_gas_stars)),1);

%figure
%[n,x] = hist(pv_dm_stars,nbins);bar(x,n./length(~isnan(pv_dm_stars)),1);
