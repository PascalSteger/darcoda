%rvir.m

figprep('virial radii distribution',...
        'log_{10}r_{vir}',...
        'log_{10}p(r_{vir})');
hold on;

histlog(hp_rvir,15);
histlog(ahf_rvir,15,'b');

legend('hp','','ahf','');