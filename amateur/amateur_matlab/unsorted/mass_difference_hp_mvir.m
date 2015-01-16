clf;

mh = hp_mvir'*1e10;
ma = ahf_mvir;

exl = log10(ma)./log10(mh)>1.007;

hold on;

scatter(log10(ma),log10(mh),'bx');
scatter(log10(ma(exl)),log10(mh(exl)),'r.');

hold off;
