%start;

hp_lxs = hp_l(1,:,1)./hp_ltot(1,:,1);
hp_lys = hp_l(2,:,1)./hp_ltot(1,:,1);
hp_lzs = hp_l(3,:,1)./hp_ltot(1,:,1);

scatter3(hp_lxs, hp_lys, hp_lzs);
%scatter3(ahf_lx, ahf_ly, ahf_lz);
axis equal;