%amautinput
%amautprepare

% reproduce fig. 10 of Croft et al. 2008

%hist(acos(cosphirl0)*180/(2*pi),50)
%hist(cosphirl0,50)
%hist(acos(cosphirl1)*180/(2*pi),50)
%hist(cosphirl0,20)
nbin=20;

[n1,xout1]=hist(acos(cosphi13)*180/(2*pi),nbin);
plot(xout1,n1,'r--')
hold on;
[n2,xout2]=hist(acos(cosphi23)*180/(2*pi),nbin);
plot(xout2,n2)
hold off;
xlabel('\phi with DM');
ylabel('p(\phi)');
legend('gas','stars');