%ahfinput
%amautinput
N=12;

hold on;
for i=1:3
    bnz=hp_lambda(:,1,i)>0 & hp_lambda(:,1,i)<5e3 & exc;
    ntot=sum(bnz);
    
    [hp_n,hp_xout]=hist(hp_lambda(bnz,1,i),N);
    figure; bar(hp_xout/10000,hp_n./ntot,1);
    if i==1
        legend('all');
    end
    if i==2
        legend('gas');
    end
    if i==3
        legend('stars');
    end;
end
hold off;