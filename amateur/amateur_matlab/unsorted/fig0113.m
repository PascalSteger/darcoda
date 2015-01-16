%ahfinput();

bnz=mvir>1e11 & mvir<1e12;
bnz=mvir>0;

XC=xc(bnz)-xc(1);
YC=yc(bnz)-yc(1);
ZC=zc(bnz)-zc(1);
LX=lx(bnz);
LY=ly(bnz);
LZ=lz(bnz);

r=sqrt(XC.*XC+YC.*YC+ZC.*ZC);
l=sqrt(LX.*LX+LY.*LY+LZ.*LZ);
cosphi=(LX.*XC+LY.*YC+LZ.*ZC)./(r.*l);

[n,xout]=hist(cosphi,20);
bar(xout,n/length(cosphi),1);
xlabel('cos \phi'); ylabel('P(m)');