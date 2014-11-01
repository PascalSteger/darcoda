%ahfinput
%amautinput

acp=cosphiauto(:,:,16);% dm stars

acp=cosphiauto(:,:,9); % gas stars

acp=cosphiauto(:,:,10);% gas dm


%acp=ahf_cpa;

%peak at 1 independent of mass
%bnz1=ahf_mvir<1e10;
%bnz1=ahf_r2<100;
%bnz1=ahf_hostno>0;
%bnz1=ahf_n_gas>0; %success. peak isolated

bnz1=ahf_n_gas>0 & exc;
bnz2=ahf_n_gas>0 & exc;

nbins=20;
subplot(1,2,1);

[n,xout]=hist(acp(bnz1),nbins);
bar(xout,n/sum(bnz1),1);
hold on;
[n,xout]=hist(abs(acp(bnz1)),nbins);
bar(xout,n/sum(bnz1),0.5,'r');
hold off;
axis([-1 1 0 0.5]);
subplot(1,2,2);

[n,xout]=hist(acp(bnz2),nbins);
bar(xout,n/sum(bnz2),1);
hold on;
[n,xout]=hist(abs(acp(bnz2)),nbins);
bar(xout,n/sum(bnz2),0.5,'r');
hold off;
axis([-1 1 0 0.5]);