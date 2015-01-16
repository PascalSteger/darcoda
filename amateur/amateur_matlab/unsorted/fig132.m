mvirpart =[];

for k=1:length(ahf_hostno)-1
    hn = ahf_hostno(k);
    if hn<=0
        continue;
    end
    mvirpart10(k) = log10(ahf_mvir(k)/ahf_mvir(hn+1));
end

ncell=20;
[n,xout] = hist(mvirpart10,ncell);
subplot(1,2,1)
bar(xout,log10(n),1);



subplot(1,2,2)
for k=1:ncell-1
    %n(k)=log10(n(k)/(log(10^xout(k+1))-log(10^xout(k))));
    %attention: use of log(10^xout) gives a too small slope:
    %-0371858x+0.39097
    %n(k)=log10(n(k)/(xout(k+1)-xout(k)));
    %-0.74201x+0.68901
    n(k)=log10(n(k)/(log(abs(xout(k)))-log(abs(xout(k+1)))));
    %-0.97648x+0.43731
end

bnz=xout>-3.5 & xout<xout(length(xout));
xout=xout(bnz);
n=n(bnz);
bar(xout,n,1);
xlabel('log_{10}(m/M)');
ylabel('log_{10}(dn/d ln(m/M))');

x=[xout', n'];
save '/data/achtland1/psteger/amd/halos/vis/report/fig132.dat' x '-ASCII';