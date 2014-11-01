% 114_hp_lambdaprime_gas_dm_stars.m: 
% lambda' parameter distribution from hp
% for dm and gas and stars

nbins = 20;
cutat = 1e100;

%% dm only
figprep;
hold on;
m_hp_lpr = m_hp_lambda(m_exc_2,1,2);
exl = m_hp_lpr > 0 & m_hp_lpr < cutat;
lpl = log10(m_hp_lpr(exl));
[n,x]=hist(lpl,nbins);
plot(x,log10(n/sum(exl)),'b');
errorbar(x,log10(n/sum(exl)),log10(n)./sqrt(n),'b.');

m_hp_lpr = m_hp_lambda(m_exc_3,1,3);
exl = m_hp_lpr > 0 & m_hp_lpr < cutat;
lpl = log10(m_hp_lpr(exl));
[n,x]=hist(lpl,nbins);
plot(x,log10(n/sum(exl)),'r');
errorbar(x,log10(n/sum(exl)),log10(n)./sqrt(n),'r.');

m_hp_lpr = m_hp_lambda(m_exc_4,1,4);
exl = m_hp_lpr > 0 & m_hp_lpr < cutat;
lpl = log10(m_hp_lpr(exl));
[n,x]=hist(lpl,nbins);
plot(x,log10(n/sum(exl)),'g');
errorbar(x,log10(n/sum(exl)),log10(n)./sqrt(n),'g.');

m_hp_lpr = m_ahf_lambda;
exl = m_hp_lpr > 0 & m_hp_lpr < cutat;
lpl = log10(m_hp_lpr(exl));
[n,x]=hist(lpl,nbins);
plot(x,log10(n/sum(exl)),'y');
errorbar(x,log10(n/sum(exl)),log10(n)./sqrt(n),'yo');

xlabel('log_{10} \lambda'); ylabel('log_{10} p(\lambda)');
%legend('dm','dm fit','gas','gas fit','stars','stars fit');
legend('dm','','gas','','stars','');
axis([-2.5,2,-2.5,0]);
hold off;

%% gas only
m_hp_lpr = m_hp_lambda(:,1,3);

exl = 0 < m_hp_lpr < cutat;

minnum = 5;
[n2, xout2]=hist(log10(m_hp_lpr(exl)),nbins);
xout2 = xout2(n2>minnum); n2=n2(n2>minnum);
errorbar(xout2,log10(n2/sum(exl)),log10(n2)./sqrt(n2),'r.');

xdata = xout2(xout2>-1.5 & xout2<0)'; ydata = log10(n2(xout2>-1.5&xout2<0)/sum(exl))'; % must be column vectors for fit input
xdata = xdata(ydata>-Inf); ydata = ydata(ydata>-Inf);
f2 = fit(xdata, ydata, fittype({'x^2','x','1'},'coefficients',{'b2','b1','b0'})); 
plot(xout2,f2(xout2), 'r-');

%% stars only
m_hp_lpr = m_hp_lambda(:,1,4);

% all halos
exl = 0 < m_hp_lpr < cutat;

minnum = 5;
[n3, xout3]=hist(log10(m_hp_lpr(exl)),nbins);
xout3=xout3(n3>minnum); n3=n3(n3>minnum);
errorbar(xout3,log10(n3/sum(exl)),log10(n3)./sqrt(n3),'g.');

xdata = xout3'; ydata = log10(n3/sum(exl))'; % must be column vectors for fit input
xdata = xdata(ydata>-Inf); ydata = ydata(ydata>-Inf);
f3 = fit(xdata, ydata, fittype({'x^2','x','1'},'coefficients',{'c2','c1','c0'})); 
plot(xout3,f3(xout3), 'g-');

%% output for SM parsing later on
%x=[xout1', n1', xout2', n2', xout3', n3'];
%save '/data/achtland1/psteger/amd/halos/vis/report/fig114.dat' x '-ASCII';
