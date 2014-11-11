figprep('distance from cluster center: hp vs ahf',...
        'd_{AHF} [h^{-1} Mpc]',...
        'd_{HP} [h^{-1} kpc]');
axis([0 300 0 3*10^5]);
set(gca,'PlotBoxAspectRatio',[1 1 1]);
plot([0:300],1000*[0:300],'g');
scatter(sqrt(ahf_pcx.*ahf_pcx+ahf_pcy.*ahf_pcy+ahf_pcz.*ahf_pcz),...
        norma(hp_xcm),...
        25*ones(size(hp_xcm(:,1))),...
        log10(hp_mtot),...
        '.')