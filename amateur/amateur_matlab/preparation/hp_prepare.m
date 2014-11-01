% hp_prepare.m: Generate the basic quantities used for an analysis of the halo properties

% array for host (or, if not present, main) halo
hostno=init_N1;
for k=1:length(ahf_hostno)
   temp=ahf_hostno(k);
   if(temp<0)        % no real host detected by merger tree
        temp=0;
   end
   hostno(k)=temp+1; % +1 since ahf_hostno starts counting with 0
end

% distance to biggest (first) halo
hp_dbig = init_N3;
hp_dbig(:,1) = hp_xcm(:,1)-hp_xcm(1,1);
hp_dbig(:,2) = hp_xcm(:,2)-hp_xcm(1,2);
hp_dbig(:,3) = hp_xcm(:,3)-hp_xcm(1,3);
hp_dbigtot   = normN3(hp_dbig);

% distance to host halo (or biggest, if none)
hp_dhost = init_N3;
hp_dhostvir = init_N3;
for k=1:length(hp_dbig)
        hp_dhost(k,:) = hp_xcm(k,:)-hp_xcm(hostno(k),:);     
        hp_dhostvir(k,:) = hp_dhost(k,:)/hp_rvir(hostno(k));
end
hp_dhosttot   = normN3(hp_dhost);
hp_dhostvirtot= normN3(hp_dhostvir);

hp_jtot = normN35(hp_j);

% arrays for basic correlations between vectors
hp_cosauto   = [];
hp_cosrj     = [];
hp_cosrea    = []; hp_cosreb    = []; hp_cosrec    = [];

% loop over all different particle types
for k=1:5
    % for autocorrelation between different halo types
    for j=1:5
        hp_cosauto= cat(3,hp_cosauto, cosN3(hp_j(:,:,k), hp_j(:,:,j)));
    end
    hp_cosrj  = cat(3, hp_cosrj,  cosN3(hp_dhost, hp_j(:,:,k)));

    hp_cosrea = cat(3, hp_cosrea, cosN3(hp_dhost, hp_ea(:,:,k)));
    hp_cosreb = cat(3, hp_cosreb, cosN3(hp_dhost, hp_eb(:,:,k)));
    hp_cosrec = cat(3, hp_cosrec, cosN3(hp_dhost, hp_ec(:,:,k)));
end
hp_cosjea = cosN35(hp_j, hp_ea);

% correlations between hp position vector and ahf angular momentum
ahf_cosrj    = cosN3(hp_dhost,cat(2,ahf_jx,ahf_jy,ahf_jz));
