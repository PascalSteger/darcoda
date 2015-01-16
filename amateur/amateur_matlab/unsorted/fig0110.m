ahfinput();

bnz=mvir_gas>0;
mvirbnz=mvir(bnz);
mvirgasbnz=mvir_gas(bnz);
fracgasvir=mvirgasbnz./mvirbnz;

[n,xout]=hist(log10(fracgasvir),50);
bar(xout,n/length(fracgasvir),1);
xlabel('log_{10}(m_{gas}/m_{vir} [msun])'); ylabel('P(m)');
