%jcorr_dm_stars_asfct_r.m

nbin = 8;
figprep('correlation for cos(j_{dm},j_{star}) as function of virial mass',...
        'log_{10}(M_{vir}/h/10^{10}M_{Sun})',...
        'cos(j_{dm},j_{star})');

exl = m_exc;
exl = exl & ~isnan(m_hp_cosauto(:,1,9));
%exl = exl & m_hp_mtot < 1;
exl = exl & m_hp_dhostvirtot < 5;

%r = m_hp_dhostvirtot(exl);
r = log10(m_hp_mvir(exl));
c = m_hp_cosauto(exl,1,9);
scatter(r,c,'b.');

p = cumsum(100/nbin*ones(nbin,1)'); %percent values
eps = 1e-10;
y = [min(r) prctile(r,p-eps)];
[n,i] = histc(r,y);
n=n(1:length(n)-1)';
mc =[]; vc =[];
for k=1:nbin
        mc(k) = mean(c(i==k));
        vc(k) = var(c(i==k));
end
errorbar(cummean(y),mc,sqrt(vc)./sqrt(n),'r.')

% smallest nonzero virial radius
%rmin = min(r(r>min(r)));
%ls = logspace(,log10(max(r)),nbin);
% count elements in intervals
% plot
%plot(ls)