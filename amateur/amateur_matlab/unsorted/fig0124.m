%amautinput
%amautprepare

% reproduction of fig. 4 Pereira et al. 2008
N=10;
nbins=10*ones(N);
rmax=10; %in r_{vir} of host halo
r=rtot; r2=r; mvir2=mvir;
acp=abs(cosphirl0(r<rmax));
%acp=abs(cosphirv0(rtot<rmax));

r2=r(r<rmax); mvir2=mvir(r<rmax);
for j=1:N
    %[nbin,ibin]=histc(rtot,linspace(0,max(rtot),nbins(j)+1));
    [nbin,ibin]=histc(r2,linspace(0,rmax,nbins(j)+1));
    rmed = []; phimed = []; phierru = []; phierrl = [];
    for i=1:nbins
        %bnzj = npart1'>200*(j-1) & npart1'<200*j & ibin==i & ~isnan(acp);
        %bnzj = npart1'>2000*(j-1) & ibin==i & ~isnan(acp) & r2<1e5;
        mres=5e10;
        bnzj = mvir2'>(j-1)*mres & ibin==i & ~isnan(acp);
        rmed(i)=-1; phimed(i)=-1;
        rmed(i) = median(r2(bnzj));
        phimed(i) = median(acp(bnzj));
        phierru(i) = prctile(acp(bnzj),84)-phimed(i);
        phierrl(i) = phimed(i) - prctile(acp(bnzj),16);
    end
    subplot(2,N/2,j)
    errorbar(rmed,phimed,phierru,phierrl,'Marker','hexagram','LineStyle','-','Color',[1 0 0]);
    axis([0,rmax,0,1]);
    hold on; plot(0.5*ones(rmax,1),'--'); hold off;

    xlabel('r / r_{vir}'); ylabel(sprintf('<|cos(r,L)|> @ m_{vir}> %d mres',(j-1)));
end