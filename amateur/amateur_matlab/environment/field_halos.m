clf


% split according to halo type

mth = [];
for k=1:length(ahf_hostno);
    if ahf_hostno(k) < 0
        mth(k) = ahf_hostno(k);
        continue;
    end
    if ahf_hostno(ahf_hostno(k)+1) < 0
        mth(k) = ahf_hostno(k);
    else
        mth(k) = -3;
    end
end

exl2 = exc_1 & ahf_n_gas > 0 & mth == -2;
exl1 = exc_1 & ahf_n_gas > 0 & mth == -1;

exl1 = exc_1 & ahf_n_gas > 0 & mth <   0;
exl0 = exc_1 & ahf_n_gas > 0 & mth >=  0;

hold all
scatter3(pcx(exl2),pcy(exl2),pcz(exl2),'r')
scatter3(pcx(exl1),pcy(exl1),pcz(exl1),'g')
scatter3(pcx(exl0),pcy(exl0),pcz(exl0),'b')
