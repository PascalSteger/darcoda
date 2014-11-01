% checking that positions of COM do not differ a lot from AHF to HP.
ha_dx = 1e3*ahf_xc'-hp0cm(1,:);
ha_dy = 1e3*ahf_yc'-hp0cm(2,:);
ha_dz = 1e3*ahf_zc'-hp0cm(3,:);

scatter3(ha_dx,ha_dy, ha_dz);

% get difference in radius from host halo
ha_rtot = sqrt(ha_dx.*ha_dx + ha_dy.*ha_dy + ha_dz.*ha_dz);

hist(log10(ha_rtot),50)

% find cosphi from ahf
ahf_dx = ahf_xc - ahf_xc(1);
ahf_dy = ahf_yc - ahf_yc(1);
ahf_dz = ahf_zc - ahf_zc(1);
ahf_dtot = sqrt(ahf_dx.*ahf_dx + ahf_dy.*ahf_dy + ahf_dz.*ahf_dz);
ahf_ltot = sqrt(ahf_lx.*ahf_lx + ahf_ly.*ahf_ly + ahf_lz.*ahf_lz);

ahf_cosphi_l = abs((ahf_lx.*ahf_dx+ahf_ly.*ahf_dy+ahf_lz.*ahf_dz)./(ahf_dtot.*ahf_ltot));
acp = ahf_cosphi_l;
hist(acp)

nbins = 10; rmax = floor(10*ahf_rvir(1)/1e3);
[nbin,ibin] = histc(ahf_dtot,linspace(0,rmax,nbins+1));

rmed1 = []; phimed1 = []; phierru1 = []; phierrl1 = [];
rmed2 = []; phimed2 = []; phierru2 = []; phierrl2 = [];
for i=1:nbins
    bnzj1 = ibin==i & ~isnan(ahf_cosphi);
    sm1 = sum(bnzj1);
    rmed1(i)=-1; phimed1(i)=-1;
    rmed1(i) = median(ahf_dtot(bnzj1));
    phimed1(i) = median(acp(bnzj1));
    phierru1(i) = (prctile(acp(bnzj1),84)-phimed1(i))/sqrt(sm1);
    phierrl1(i) = (phimed1(i) - prctile(acp(bnzj1),16))/sqrt(sm1);
end
errorbar(rmed1,phimed1,phierru1,phierrl1,'Marker','hexagram','LineStyle','-','Color',[1 0 0]);
axis([0,rmax,0,1]);
hold on; plot(0.5*ones(rmax,1),'--'); hold off;
xlabel('r / r_{vir}'); ylabel(sprintf('<|cos(r,L)|>'));
