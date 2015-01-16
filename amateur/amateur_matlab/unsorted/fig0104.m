ahfinput();

scale=35;
border=10/scale;
bnz= mvir>1e10 & mvir<1e11 & lambda>0 & lambda<border;
lambda=lambda(bnz);

[n,xout]=hist(lambda(lambda<border),scale);
bar(xout,n/length(lambda(lambda<border)),1); xlabel('\lambda'); ylabel('P(\lambda)');
hold on;
%histfit(log10(lambda(lambda<border)),20);
%return;

%to 5% certainty
[parmhat,parmci] = lognfit(lambda,0.05);
exp(parmhat(1))
exp(parmci(:,1))
x=[0:border/100:border];
plot(x,lognpdf(x,parmhat(1),parmhat(2))/(4*scale));
hold off;