figprep('physical position','x','y');
hold on;
S = ones(size(ahf_pcx));
sc = 1000; % difference between AHF (in Mpc) and hp (in kpc)
scatter3(ahf_pcx*sc,ahf_pcy*sc,ahf_pcz*sc,30*S,'r')
scatter3(hp_xcmx,hp_xcmy,hp_xcmz,10*S,'b')
