% fig140_gasfractions: histograms of gas properties
clf


subplot(1,2,1);set(gca,'FontSize',15);
set(gca,'FontSize',15);
[n1,x1]=hist(log10(snap_uint),500);
stairs(x1,n1/length(snap_uint),'b');

hold on;

[n2,x2]=hist(log10(snap_temp),500);
stairs(x2,n2/length(snap_temp),'r');
xlabel('log_{10}u_{gas},log_{10}T_{gas}');
ylabel('frequency');

hold off;

%%
%loglog(snap_temp,snap_uint./snap_temp,'.')
%xlabel('log_{10}T'); ylabel('log_{10}n');

%% gas fraction of all halos as function of ...

% total mass of haloes with a min. amount of gas as given in exclude.m
subplot(1,2,2);set(gca,'FontSize',15);
loglog(m_hp_mass(logical(m_exc_2),1,1),m_hp_mass(logical(m_exc_2),1,2)./m_hp_mass(logical(m_exc_2),1,1),'.'); xlabel('total halo mass');ylabel('gas fraction');

% radius
%loglog(r(logical(m_exc_2)),m_hp_mass(logical(m_exc_2),1,2)./m_hp_mass(logical(m_exc_2),1,1),'.'); xlabel('r_{vir}');ylabel('gas fraction');