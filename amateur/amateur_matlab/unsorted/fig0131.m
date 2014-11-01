%amautinput
%amautprepare
%tidalfieldinput

nbins=10;
ewmax=max(ew1);
rmax=3;
r=rtot; r2=r(r<rmax);

cosphi1rt=(cmxold.*t1x+cmyold.*t1y+cmzold.*t1z)./(t1tot.*rtotold);
cosphi2rt=(cmxold.*t2x+cmyold.*t2y+cmzold.*t2z)./(t2tot.*rtotold);
cosphi3rt=(cmxold.*t3x+cmyold.*t3y+cmzold.*t3z)./(t3tot.*rtotold);

cosphi01lt=(hp0lx.*t1x+hp0ly.*t1y+hp0lz.*t1z)./(t1tot.*hp0ltot);
cosphi02lt=(hp0lx.*t2x+hp0ly.*t2y+hp0lz.*t2z)./(t2tot.*hp0ltot);
cosphi03lt=(hp0lx.*t3x+hp0ly.*t3y+hp0lz.*t3z)./(t3tot.*hp0ltot);

cosphi11lt=(hp1lx.*t1x+hp1ly.*t1y+hp1lz.*t1z)./(t1tot.*hp1ltot);
cosphi12lt=(hp1lx.*t2x+hp1ly.*t2y+hp1lz.*t2z)./(t2tot.*hp1ltot);
cosphi13lt=(hp1lx.*t3x+hp1ly.*t3y+hp1lz.*t3z)./(t3tot.*hp1ltot);

cosphi21lt=(hp2lx.*t1x+hp2ly.*t1y+hp2lz.*t1z)./(t1tot.*hp2ltot);
cosphi22lt=(hp2lx.*t2x+hp2ly.*t2y+hp2lz.*t2z)./(t2tot.*hp2ltot);
cosphi23lt=(hp2lx.*t3x+hp2ly.*t3y+hp2lz.*t3z)./(t3tot.*hp2ltot);

cosphi31lt=(hp3lx.*t1x+hp3ly.*t1y+hp3lz.*t1z)./(t1tot.*hp3ltot);
cosphi32lt=(hp3lx.*t2x+hp3ly.*t2y+hp3lz.*t2z)./(t2tot.*hp3ltot);
cosphi33lt=(hp3lx.*t3x+hp3ly.*t3y+hp3lz.*t3z)./(t3tot.*hp3ltot);

acp=abs(cosphi1rt(r<rmax));

mvir2=mvir(r<rmax); mres=1e11;

%[nbin,ibin]=histc(rtot,linspace(0,max(rtot),nbins(j)+1));
[nbin,ibin]=histc(r2,prctile(r2,[linspace(0,100,nbins+1)]));
rmed1 = []; phimed1 = []; phierru1 = []; phierrl1 = [];
rmed2 = []; phimed2 = []; phierru2 = []; phierrl2 = [];
nnzj1 = []; nnzj2 = [];
for i=1:nbins    
    bnzj1 = mvir2'>10*min(mvir(npart1==min(npart1))) & mvir2'<mres & ibin==i & ~isnan(acp);
    nnzj1(i) = sum(bnzj1);
    
    bnzj2 = mvir2'>mres & ibin==i & ~isnan(acp);
    nnzj2(i) = sum(bnzj2);
    
    rmed1(i)=-1; phimed1(i)=-1;
    rmed2(i)=-1; phimed2(i)=-1;
    rmed1(i) = median(r2(bnzj1));
    rmed2(i) = median(r2(bnzj2));
    phimed1(i) = mean(acp(bnzj1));
    phimed2(i) = mean(acp(bnzj2));
    phierru1(i) = prctile(acp(bnzj1),84)-phimed1(i);
    phierru2(i) = prctile(acp(bnzj2),84)-phimed2(i);
    phierrl1(i) = phimed1(i) - prctile(acp(bnzj1),16);
    phierrl2(i) = phimed2(i) - prctile(acp(bnzj2),16);
end
subplot(1,2,1)
errorbar(rmed1,phimed1,phierru1,phierrl1,'Marker','hexagram','LineStyle','-','Color',[1 0 0]);
hold on;
scatter(r2(mvir2'>10*min(mvir(npart1==min(npart1))) & mvir2'<mres),acp(mvir2'>10*min(mvir(npart1==min(npart1))) & mvir2'<mres));
axis([0,rmax,0,1]);
 plot(0.5*ones(rmax,1),'--'); hold off;
xlabel('r/r_{vir}'); ylabel(sprintf('<|cos(r,t)|> @ m_{vir}< %d',mres));

subplot(1,2,2)
errorbar(rmed2,phimed2,phierru2,phierrl2,'Marker','hexagram','LineStyle','-','Color',[1 0 0]);
hold on;
axis([0,rmax,0,1]);
scatter(r2(mvir2'>mres),acp(mvir2'>mres));
 plot(0.5*ones(rmax,1),'--'); hold off;
xlabel('r/r_{vir}'); ylabel(sprintf('<|cos(t,t)|> @ m_{vir}> %d',mres));
