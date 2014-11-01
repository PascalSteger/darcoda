%fig0116
% subhalo mass function
%ahfinput()
%amautinput

hostids=ahf_subno;
mvirpart=[]; mvirpart10=[];
for k=1:length(ahf_mvir)
    if ahf_mvir(k)==0
        continue;
    end
    if ahf_mvir(k)>=max(ahf_mvir)
       continue;
    end
    if ahf_hostno(k)<0
        continue;
    end
    if(ahf_npart1(k)/ahf_npart1(ahf_hostno(k)+1)>0.7)
        continue;
    end

    mvirpart10(k) = log10(ahf_mvir(k)/ahf_mvir(ahf_hostno(k)+1));
    mvirpart(k) = log(ahf_mvir(k)/ahf_mvir(ahf_hostno(k)+1));
end

%bnz = mvir>0 & mvir~=max(mvir);
%mvirpart10 = log10(mvir(bnz)/max(mvir));
%mvirpart = log(mvir(bnz)/max(mvir));

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
    n(k)=log10(n(k)/(log(xout(k+1))-log(xout(k))));
    %-0.97648x+0.43731
end

bnz=xout>-3.5 & n>-Inf;
xout=xout(bnz);
n=n(bnz);
bar(xout,n,1);
xlabel('log_{10}(m/M)');
ylabel('log_{10}(dn/d ln(m/M))');

