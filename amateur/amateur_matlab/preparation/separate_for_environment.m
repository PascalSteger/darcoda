% select all halos belonging to clusters
%          lambda1,2,3>=0:           all eigenvalues positive, stable orbit

exc_cluster = (Tew(1,:) >= 0 & Tew(2,:) >= 0 & Tew(3,:) >= 0);

% belonging to filaments
%          lambda2,3>=0;    lambda1<0

exc_filament = (Tew(1,:) < 0 & Tew(2,:) >= 0 & Tew(3,:) >= 0);

% belonging to sheets
%          lambda3>=0;      lambda1,2<0

exc_sheet = (Tew(1,:) < 0 & Tew(2,:) < 0 & Tew(3,:) >= 0);

% belonging to voids
%          lambda1,2,3<0

exc_void = (Tew(1,:) < 0 & Tew(2,:) < 0 & Tew(3,:) < 0);

%exc = exc_cluster & hp_pcx;
%p1=plot3(hp_pcx(exc_cluster),hp_pcy(exc_cluster),hp_pcz(exc_cluster),'or');
%hold all;
%p2=plot3(hp_pcx(exc_filament),hp_pcy(exc_filament),hp_pcz(exc_filament),'.b');

%plot3(hp_pcx(exc_sheet),hp_pcy(exc_sheet),hp_pcz(exc_sheet),'o');


% simple example of subhalo mass function

%hold off;
%[hc,xc]=hist(log(hp_mvir(exc_cluster)),20);
%[hf,xf]=hist(log(hp_mvir(exc_filament)),20);
%[hs,xs]=hist(log(hp_mvir(exc_sheet)),20);
%[hv,xv]=hist(log(hp_mvir(exc_void)),10);

%hold all;
%stairs(xc,hc/sum(exc_cluster),'b')
%stairs(xf,hf/sum(exc_filament),'r')
%stairs(xs,hs/sum(exc_sheet),'g')
%stairs(xv,hv/sum(exc_void),'black')
