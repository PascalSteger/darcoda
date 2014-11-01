%ahfinput();

r=sqrt(xc.*xc+yc.*yc+zc.*zc);
l_gas=sqrt(lx_gas.*lx_gas+ly_gas.*ly_gas+lz_gas.*lz_gas);
cosphi=(lx_gas.*xc+ly_gas.*yc+lz_gas.*zc)./(r.*l_gas);

[n,xout]=hist(cosphi,20);
bar(xout,n/length(cosphi),1);
xlabel('cos \phi'); ylabel('P(m)');