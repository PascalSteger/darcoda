ptype=hdf5read('/data/branwyn1/psteger/AMD/snap/snapshot_008.hdf5','Data_Type');
pos=hdf5read('/data/branwyn1/psteger/AMD/snap/snapshot_008.hdf5','Data_Pos');
vel=hdf5read('/data/branwyn1/psteger/AMD/snap/snapshot_008.hdf5','Data_Vel');
mass=hdf5read('/data/branwyn1/psteger/AMD/snap/snapshot_008.hdf5','Data_Mass');

bnz=(mass>1000);

X=pos(1,:); Y=pos(2,:); Z=pos(3,:);

%scatter3(X(bnz),Y(bnz),Z(bnz));
hist(ptype,5,1);