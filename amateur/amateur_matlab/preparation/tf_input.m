% tf_input.m: Get information from TidalField

Tpos=[];
Tew=[]; Teva=[]; Tevb=[]; Tevc=[];
Tdens=[]; Tovdens=[]; Tpot=[]; Ttens=[];

Tpos=hdf5read(file_tensor,'Coordinates');
Tew=hdf5read(file_tensor,'Eigenvalues')';
Teva=hdf5read(file_tensor,'Eigenvector_1')';
Tevb=hdf5read(file_tensor,'Eigenvector_2')';
Tevc=hdf5read(file_tensor,'Eigenvector_3')';

Tdens=hdf5read(file_tensor,'Density');
Tovdens=hdf5read(file_tensor,'Overdensity');
Tpot=hdf5read(file_tensor,'Potential');
Ttens=hdf5read(file_tensor,'TensorField');

% eigenwerte in Tew(123,:)
tt=cat(3,Teva,Tevb,Tevc); %xyz | halo number | abc(ev number)
ttot=[];
% should be normalized to 1, already done by TidalField. To check this:
for k=1:3
    ttot=cat(3,ttot, normN3(tt(:,:,k)));
end
