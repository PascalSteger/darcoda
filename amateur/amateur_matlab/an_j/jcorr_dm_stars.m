%jcorr_dm_stars.m

nbin = 20;
figprep('correlation for cos(j_{dm},j_{star})',...
        '\xi:=cos(j_{dm},j_{star})',...
        'p(\xi)');

exl = m_exc;
[n,x] = hist(m_cosauto(exl,1,9),nbin);
y = n/sum(n);
plot(x,y,'r');
errorbar(x,y,y./sqrt(n),'r.');
flat = [-1:0.01:1];
plot(flat,ones(size(flat))/nbin)