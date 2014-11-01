ahfinput();

[n,xout]=hist(log10(mvir),50);
bar(xout,n/length(mvir),1);
xlabel('log_{10}(m_{vir}/m_{Sun})'); ylabel('P(m)');
