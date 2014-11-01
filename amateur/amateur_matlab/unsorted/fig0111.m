ahfinput();

r=sqrt(xc.*xc+yc.*yc+zc.*zc);
l=sqrt(lx.*lx+ly.*ly+lz.*lz);
cosphi=(lx.*xc+ly.*yc+lz.*zc)./(r.*l);

[n,xout]=hist(cosphi,20);
bar(xout,n/length(cosphi),1);
xlabel('cos \phi'); ylabel('P(m)');