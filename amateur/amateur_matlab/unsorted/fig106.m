% aim: reproduce fig.2 of Pereira, Bryan and Gill 2008
% here hp is used
start;

% major axis is given by vector from biggest eigenvalue,
% use 

angle_r_em = cosine(hp_dxc,hp_dyc,hp_dzc,hp_ea(1,:,4),hp_ea(2,:,4),hp_ea(3,:,4));
%angle_r_em = cosine(hp_dxc,hp_dyc,hp_dzc,hp_ea(1,:,4),hp_ea(2,:,4),hp_ea(3,:,4));
%angle_r_em = cosine(ahf_dxc,ahf_dyc,ahf_dzc,hp_em(1,:,1),hp_em(2,:,1),hp_em(3,:,1));
% for all constituents, take 1,:,1/6
%angle_r_em = cosphireaold(1,:,4);

[n,x] = hist(abs(angle_r_em),15);

bar(x,n/sum(~isnan(angle_r_em)),1);
xlabel('cos(r, I)'); ylabel('fraction');