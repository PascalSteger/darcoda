ahfinput();

scale=20;
border=10/scale;
bnz= mvir>1e11 & mvir<1e12 & lambda_gas>0 & lambda_gas<border;
lambda_gas=lambda_gas(bnz);

[n,xout]=hist(lambda_gas(lambda_gas<border),scale);
bar(xout,n/length(lambda_gas(lambda_gas<border)),1); xlabel('\lambda_{gas}'); ylabel('P(\lambda_{gas})');
hold on;
%histfit(log10(lambda_gas(lambda_gas<border)),20);
%return;

%to 5% certainty
[parmhat,parmci] = lognfit(lambda_gas,0.05);
exp(parmhat(1))
exp(parmci(:,1))
x=[0:border/100:border];
plot(x,lognpdf(x,parmhat(1),parmhat(2))/(4*scale));
hold off;