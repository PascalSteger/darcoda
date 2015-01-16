%ahfinput
%amautinput

acp=abs(cosphiauto(:,:,10));
bnzacp=acp>0.95;

%peak at 1 independent of mass
%bnz1=ahf_mvir<1e10;
%bnz1=ahf_r2<100;
%bnz1=ahf_hostno>0;
%bnz1=ahf_n_gas>0; %success. peak isolated
bnz1=ahf_mvir_gas>0;
bnz2=~bnz1;

nbins=40;
subplot(1,2,1);
axis tight;
[n,xout]=hist(acp(bnz1),nbins);
bar(xout,n/sum(bnz1));

subplot(1,2,2);
axis tight;
[n,xout]=hist(acp(bnz2),nbins);
bar(xout,n/sum(bnz2));