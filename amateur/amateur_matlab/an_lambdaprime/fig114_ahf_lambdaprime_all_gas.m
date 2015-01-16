% 114_lambdaprime.m: 
% lambda' parameter distribution from ahf and hp

nbins = 10;
cutat = 1e100;

%% ahf only
clf;hold on;
lambda = m_ahf_lambda;
exl = lambda > 0 & lambda < cutat;
lpl = log10(lambda(exl));
[n,x]=hist(lpl,nbins);
%scatter(x,log10(n/sum(exl)));
errorbar(x,log10(n/sum(exl)),log10(n)./sqrt(n),'b.');

%% hp only
lambda = m_hp_lambda(:,1,1);
exl = lambda > 0 & lambda < cutat;
lpl = log10(lambda(exl));
[n,x]=hist(lpl,nbins);
%scatter(x,log10(n/sum(exl)));
errorbar(x,log10(n/sum(exl)),log10(n)./sqrt(n),'b.');


%% output for SM parsing later on
%x=[xout1', n1', xout2', n2', xout3', n3'];
%save '/data/achtland1/psteger/amd/halos/vis/report/fig114.dat' x '-ASCII';
