%0116_subhmfct_linreg_lowmasscut.m: subhalo mass function
figprep('subhalo mass function',...
        'log_{10}m/M',...
        'log_{10}(dn/d ln(m/M))');

% exclude all subhaloes that are not primary
mth = [];
for k=1:length(m_ahf_hostno);
    if m_ahf_hostno(k) < 0
        mth(k) = m_ahf_hostno(k);
        continue;
    end
    if m_ahf_hostno(m_ahf_hostno(k)+1) < 0
        mth(k) = m_ahf_hostno(k);
    else
        mth(k) = -3;
    end
end

m_hostids = m_ahf_subno;
m_mvirpart = []; m_mvirpart10 = [];
for k=1:length(m_ahf_mvir)
    % if we have too little particles, let's continue
    if m_exc_1(k) == 0
        m_mvirpart(k) = NaN;
        m_mvirpart10(k) = NaN;
        continue;
    end

    % exclude haloes with no real host as well
    if m_ahf_hostno(k) < 0
        m_mvirpart(k) = NaN;
        m_mvirpart10(k) = NaN;
        continue;
    end
    % exclude haloes that share too many particles with their host
    % i.e., exclude "fake" subhaloes
    %if(m_ahf_npart1(k)/m_ahf_npart1(m_ahf_hostno(k)+1)>0.7)
    %    continue;
    %end
    % calculate the SHMF using ahf...
    %m_mvirpart(k) = log(m_ahf_mvir(k)/m_ahf_mvir(m_ahf_hostno(k)+1));
    m_mvirpart10(k) = log10(m_ahf_mvir(k)/m_ahf_mvir(m_ahf_offset(k)+1));
    % or hp properties
    %m_mvirpart10(k) = log10(m_hp_mtot(k)/m_hp_mtot(1));
    m_mvirpart10(k) = log10(m_hp_mtot(k)/m_hp_mtot(m_ahf_offset(k)+1));
    %m_mvirpart10(k) = log10(m_hp_mtot(k)/m_hp_mtot(m_ahf_hostno(k)+1));
    
end

% if we want to plot all range
ncell = 20;
[n,xout] = hist(m_mvirpart10,ncell);
%subplot(1,2,1);set(gca,'FontSize',15);
%bar(xout(1:length(xout)),log10(n(1:length(n))),1);
%xlabel('log_{10}(m/M)');
%ylabel('log_{10}(dn/d ln(m/M))');

% cut low mass end (numeric effect: +1 particle gives too big a change)

for k=1:ncell-1
    %n(k)=log10(n(k)/(log(10^xout(k+1))-log(10^xout(k))));
    %attention: use of log(10^xout) gives a slope that is too small:
    %n(k)=log10(n(k)/(xout(k+1)-xout(k)));
    n(k)=log10(n(k)/(log(xout(k+1))-log(xout(k))));
end

bnz=xout>-4 & n>-Inf;
xout=xout(bnz); n=n(bnz);

bar(xout(1:length(xout)),n(1:length(n)),1);
%bar(xout,n,1);
