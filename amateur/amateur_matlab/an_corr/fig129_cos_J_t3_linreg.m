% 129_cos_J_t3_linreg.m: correlation between angular momentum and third eigenvector of tidal field

acp = abs(m_hp_cosphijt(3,:,1));   %abc,:,1..4(6) we can choose m_(ahf)_cosphilt
%acp = abs(m_ahf_cosphitea(1,:,1)); %1,:,abc(of Teva/b/c)

ftemp = prctile(hp_mvir(exc),[100/3,200/3]);
excf = logical(exc) & hp_mvir' < ftemp(1);
excm = logical(exc) & hp_mvir' >= ftemp(1) & hp_mvir' <= ftemp(2);
excl = logical(exc) & hp_mvir' > ftemp(2);
acpf = acp(excf);
acpm = acp(excm);
acpl = acp(excl);

nbins = 10;
%[n,x] = hist(acpf,nbins);
%[n,x] = hist(acpl,nbins);
[n,x] = hist(acpl,nbins);

pn = double(n)/sum(n);
x1 = x; n1 = pn;
errorbar(x,pn,pn./sqrt(n));hold on;
%lsline;
plot([min(acp),max(acp)],[1/nbins 1/nbins],'.-.');

[p,ErrorEst] = polyfit(x,pn,1);
[pop_fit,delta] = polyval(p,x,ErrorEst);
% Plot the data, the fit, and the confidence bounds
plot(x,pn,'o',...
     x,pop_fit,'g-',...
     x,pop_fit+delta,'r--',...
     x,pop_fit-delta,'r--'); 
% Annotate the plot
% TBD: change according to first choice in line 4!
xlabel('cos(J,t_1)');
xlabel('cos(e_a,t_1)');
ylabel('percentage');
grid on;
hold off;

% output to file for later SuperMongo treatment
%x2 = x; n2 = pop_fit; e1 = delta;
%xx=[x1', n1', x2', n2', delta'];
%save '/data/achtland1/psteger/amd/halos/vis/report/fig129.dat' xx '-ASCII';

%alignment between l_gas, l_stars
%hist(abs(cosphiauto(1,exc,9)))
