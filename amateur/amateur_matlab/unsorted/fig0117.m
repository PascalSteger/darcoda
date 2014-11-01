vxc=vxc-vxc(1);
vyc=vyc-vyc(1);
vzc=vzc-vzc(1);
dx=xc-xc(1);
dy=yc-yc(1);
dz=zc-zc(1);

d=sqrt(dx.*dx+dy.*dy+dz.*dz);
v=sqrt(vxc.*vxc+vyc.*vyc+vzc.*vzc);

bnz=d>0 & v>0 & mvir>0 & log10(mvir)>-1 & log10(d)>-1;
p=plot(log10(mvir(bnz)),log10(v(bnz)),'r.');
plot(log10(d(bnz)),log10(v(bnz)),'rx');
xlabel('log_{10}(d), d>0.1Mpc');
ylabel('log_{10}(v)');
h=lsline;
set(h,'color','b');
