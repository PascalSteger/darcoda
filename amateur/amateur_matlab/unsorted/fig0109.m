ahfinput();

bnz=mvir_gas>0;
mvirgasbnz=mvir_gas(bnz);

[n,xout]=hist(log10(mvirgasbnz),50);
bar(xout,n/length(mvirgasbnz),1);
xlabel('log_{10}(m_{vir,gas}/m_{Sun})'); ylabel('P(m)');
