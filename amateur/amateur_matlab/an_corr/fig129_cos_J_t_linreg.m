% 129_cos_J_t3_linreg.m: correlation between angular momentum and third eigenvector of tidal field
clf;
set(gca,'FontSize',15);

m_acp = (m_ahf_cosjt(1,:,1));   %(Tabc,:,halo type)
%m_acp = (m_ahf_costea(1,:,1));  %(1,:,Tevabc)
% %m_acp = acp(logical(m_exc));

ftemp = prctile(m_hp_mvir(logical(m_exc_1)),[100/3,200/3]);
m_excf = logical(m_exc_1) & m_hp_mvir < ftemp(1);
m_excm = logical(m_exc_1) & m_hp_mvir >= ftemp(1) & m_hp_mvir <= ftemp(2);
m_excl = logical(m_exc_1) & m_hp_mvir > ftemp(2);
m_acpf = m_acp(m_excf); %first 1/3
m_acpm = m_acp(m_excm); %middle
m_acpl = m_acp(m_excl); %last

%KNOB
nbins = 9;
[n,x] = hist(m_acpf,nbins);
[n,x] = hist(m_acpm,nbins);
[n,x] = hist(m_acpl,nbins);
[n,x] = hist(m_acp,nbins);

pn = double(n)/sum(n);
x1 = x; n1 = pn;
errorbar(x,pn,pn./sqrt(n));
xlabel('\phi (t_1, J)'); ylabel('frequency');
%xlabel('\phi (t_1,e_1)'); ylabel('frequency');
%hold on;
%plot([min(acp),max(acp)],[1/nbins 1/nbins],'.-.');

%% fit using exp on [-1,1], using linear composition of several functions
% xdata = x'; ydata = pn'; % must be column vectors for fit input
% f = fit(xdata, ydata, fittype({'exp(x)','1'},'coefficients',{'a1','a2'})); 
% plot(f, 'r');

%% polyfit
%[p,ErrorEst] = polyfit(x,pn,2);
%[pop_fit,delta] = polyval(p,x,ErrorEst);
% Plot the data, the fit, and the confidence bounds
%plot(x,pn,'o',...
%     x,pop_fit,'g-',...
%     x,pop_fit+delta,'r--',...
%     x,pop_fit-delta,'r--'); 
% Annotate the plot
% TBD: change according to first choice in line 4!
%xlabel('cos(J,t_1)');
%xlabel('cos(e_a,t_1)');
%ylabel('percentage');
%grid on;
%hold off;
%% end polyfit

%% output to file for later SuperMongo treatment
%x2 = x; n2 = pop_fit; e1 = delta;
%xx=[x1', n1', x2', n2', delta'];
%save '/data/achtland1/psteger/amd/halos/vis/report/fig129.dat' xx '-ASCII';

%alignment between l_gas, l_stars
%hist(abs(cosauto(1,exc,9)))
