clf;
hold on;
grid on;

nbins = 40;

[n,xout]=hist(log10(m_ahf_sigv),nbins);
plot(xout,log10(n),'b.-')

[n,xout]=hist(log10(m_hp_sigma(:,1,1)),nbins);
plot(xout,log10(n),'r.-')

legend('ahf','hp')
xlabel('log_{10}\sigma');
ylabel('count');