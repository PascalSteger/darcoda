% lambda' parameter distribution from hp for dm and gas
%start;

%figure
nbins = 20;
cutat = 1e100;

%DM only
hp_lpr = m_hp_lambda(:,1,4);
%hp_lpr = m_ahf_lambda;
exl = hp_lpr < cutat & hp_lpr>0 & m_ahf_hostno'<0; sum(exl)

[n1, xout1]=hist(log10(hp_lpr(exl)),nbins);
stairs(xout1,log10(n1/sum(exl)));

hold on;

%gas only
hp_lpr = m_hp_lambda(:,1,2);
exl = hp_lpr < cutat & hp_lpr > 0; sum(exl)

[n2, xout2]=hist(log10(hp_lpr(exl)),nbins);
stairs(xout2,log10(n2/sum(exl)),'r');

%stars only
hp_lpr = m_hp_lambda(:,1,3);
exl = hp_lpr < cutat & hp_lpr>0; sum(exl)

[n3, xout3]=hist(log10(hp_lpr(exl)),nbins);
stairs(xout3,log10(n3/sum(exl)),'g');


hold off;
xlabel('log_{10} \lambda'); ylabel('log_{10} p(\lambda)');
legend('DM','gas','stars');

%x=[xout1', n1', xout2', n2', xout3', n3'];
%save '/data/achtland1/psteger/amd/halos/vis/report/fig114.dat' x '-ASCII';