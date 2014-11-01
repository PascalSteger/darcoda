data_id=hdf5read('/data/achtland1/psteger/amd/output/tarkin_25174_CSF_45x_snap_302/snapshot.hdf5','Data_ID');
data_pos=hdf5read('/data/achtland1/psteger/amd/output/tarkin_25174_CSF_45x_snap_302/snapshot.hdf5','Data_Pos');
data_vel=hdf5read('/data/achtland1/psteger/amd/output/tarkin_25174_CSF_45x_snap_302/snapshot.hdf5','Data_Vel');
data_mass=hdf5read('/data/achtland1/psteger/amd/output/tarkin_25174_CSF_45x_snap_302/snapshot.hdf5','Data_Mass');

bnz=(data_mass>1000);

X=data_pos(1,:); Y=data_pos(2,:); Z=data_pos(3,:);

%scatter3(X(bnz),Y(bnz),Z(bnz));
hist(data_id,5,1);
