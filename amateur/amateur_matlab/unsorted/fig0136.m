%ahfinput
%amautinput (prepare not required)

bnz=ahf_lambda>0 & ahf_lambda<0.2;
bnz_gas=ahf_lambda_gas>0 & ahf_lambda_gas<0.2;

ahf_lambdaloc=ahf_lambda(bnz);
ahf_lambdagasloc=ahf_lambda_gas(bnz_gas);

N=20; ntot=length(ahf_lambdaloc); ntotgas=length(ahf_lambdagasloc);
[ahf_n,ahf_xout]=hist(ahf_lambdaloc,N);
[ahf_n_gas,ahf_xout_gas]=hist(ahf_lambdagasloc,N);

bar(ahf_xout,ahf_n./ntot,'r');
hold on;
bar(ahf_xout_gas,ahf_n_gas./ntotgas,'b');
legend('DM and gas','gas');
hold off;
