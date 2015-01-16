% lambda' parameter distribution from ahf for dm and gas
%start;
exl = ahf_lambda < 10000;

[n, xout]=hist(log10(ahf_lambda(exl)),40);
stairs(xout,log10(n/sum(exl)));

hold on;

% most lambda's for ahf are 0, exclude them (no gas at all in halo)
exl = exc & ahf_lambda_gas > 0 & ahf_lambda_gas < 0.4;
[n, xout]=hist(log10(ahf_lambda_gas(exl)),20);
stairs(xout,log10(n/sum(exl)),'r');

hold off;
xlabel('log_{10}\lambda'); ylabel('log_{10}p(\lambda)');
legend('DM','gas');