% tf_prepare.m: prepare TidalField input

% needed arrays
hp_cosjt = init_N3; hp_cosrt = []; ahf_cosjt = init_N3;
hp_costea = []; hp_costeb = []; hp_costec = [];
ahf_costea = [];ahf_costeb = [];ahf_costec = [];

% loop over all eigenvectors (1,2,3)
for j=1:3
    % loop over all halo types
    hp_cosjt(:,j) = cosN3(hp_j(:,:,1),tt(:,:,j));
    ahf_cosjt(:,j) = ahf_jx.*tt(:,1,j)+ahf_jy.*tt(:,2,j)+ahf_jz.*tt(:,3,j);
    hp_cosrt = cat(3, hp_cosrt,         cosN3(hp_dhost,tt(:,:,j)));
    ahf_costea = cat(3,ahf_costea,ahf_eax.*tt(:,1,j)+ahf_eay.*tt(:,2,j)+ahf_eaz.*tt(:,3,j));
    ahf_costeb = cat(3,ahf_costeb,ahf_ebx.*tt(:,1,j)+ahf_eby.*tt(:,2,j)+ahf_ebz.*tt(:,3,j));
    ahf_costec = cat(3,ahf_costec,ahf_ecx.*tt(:,1,j)+ahf_ecy.*tt(:,2,j)+ahf_ecz.*tt(:,3,j));
end
