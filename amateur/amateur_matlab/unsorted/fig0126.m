%fig0126
%ahfinput()

bnz=mvir>0 & mvir<max(mvir);
mvirpart10=log10(mvir(bnz)/max(mvir));
mvirpart=log(mvir(bnz)/max(mvir));

ncell=40;
[n,xout] = hist(mvirpart10,ncell);
subplot(1,2,1)
stairs(xout,log10(n));

subplot(1,2,2)
for i=1:ncell-1
    %n(i)=log10(n(i)/(log(10^xout(i+1))-log(10^xout(i))));
    %attention: use of log(10^xout) gives a too small slope:
    %-0371858x+0.39097
    %n(i)=log10(n(i)/(xout(i+1)-xout(i)));
    %-0.74201x+0.68901
    n(i)=log10(n(i)/(log(xout(i+1))-log(xout(i))));
    %-0.97648x+0.43731
end

bnz=xout>log10((max(mvir(npart1==min(npart1)))-min(mvir(npart1==min(npart1))))/mvir(1)) & n>-Inf;
xout=xout(bnz);
n=n(bnz);
stairs(xout,n);
xlabel('log_{10}(m/M)');
ylabel('log_{10}(dn/d ln(m/M))');

