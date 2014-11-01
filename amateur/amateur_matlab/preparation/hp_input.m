% hp_input: Get information from HaloProperties

% initialize all used arrays
tic

A = importdata(file_hprop);

% all particles
%hp_ptypn = cat(3,hp_ptypn,  hp_pnum-hp_ptypn(:,1,1)
%-hp_ptypn(:,1,2)-hp_ptypn(:,1,3)-hp_ptypn(:,1,4));
%hp_ahf_xc = hdf5read(file_extra,'ahf_xc');

hp_ntot   = A(:,2);
hp_mtot   = A(:,3);
hp_rtot   = A(:,4);
hp_mvir   = A(:,5);
hp_rvir   = A(:,6);
hp_xcm    = A(:,7:9);
hp_vcm    = A(:,10:12);

%specific subhalos
tmp = zeros(size(A,1),3,5);
hp_j=tmp; hp_ea=tmp; hp_eb=tmp; hp_ec=tmp; 
tmp2 = zeros(size(A,1),1,5);
hp_npart=tmp2;  hp_mpart=tmp2; 
hp_lambda=tmp2; hp_sigma=tmp2;
hp_eatot=tmp2;  hp_ebtot=tmp2;  hp_ectot=tmp2;

offset=19;
for k=0:4
    hp_npart(:,:,k+1) = A(:,19+k*offset);
    hp_mpart(:,:,k+1)  = A(:,20+k*offset);
    hp_lambda(:,:,k+1)= A(:,21+k*offset);
    hp_sigma(:,:,k+1) = A(:,22+k*offset);

    hp_j(:,:,k+1)     = A(:,(23:25)+k*offset);
    hp_eatot(:,1,k+1) = A(:,26+k*offset);
    hp_ea(:,:,k+1)    = A(:,(27:29)+k*offset);
    hp_ebtot(:,1,k+1) = A(:,30+k*offset);
    hp_eb(:,:,k+1)    = A(:,(31:33)+k*offset);
    hp_ectot(:,1,k+1) = A(:,34+k*offset);
    hp_ec(:,:,k+1)    = A(:,(35:37)+k*offset);

end
toc
