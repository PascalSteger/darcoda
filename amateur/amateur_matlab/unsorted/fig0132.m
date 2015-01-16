%ahfinput();

lambdas=hp0lambda;
bnz= hp0mass>0 & lambdas>0 & lambdas<2000 & ~isnan(lambdas);
lambdas=lambdas(bnz);

[n,xout]=hist(lambdas,20);
bar(xout,n/length(lambdas),1); xlabel('\lambda'); ylabel('P(\lambda)');

%histfit(log10(lambda2(lambda2<border)),20);
%return;