%ahfinput

plot(log10(npart2),'b')
hold on;
plot(log10(npart1),'r')
plot(log10(npart2)-log10(npart1)-2,'g')
hold off;
xlabel('halo number');
ylabel('log_{10}(count)');
legend('gas','dark matter','difference');

%hist(log10(npart1),50)
%hold on;
%hist(log10(npart2),50)
%hold off;