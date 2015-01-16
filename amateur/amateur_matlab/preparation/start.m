% prepare all data for later analysis
% reads in files from ahf, hp and tf runs.
% ahf_*: AHFstep output
% hp_* : HaloProperties output
% tf_* : TidalField output

clear;
tic;

% define the numbers used for different particle and halo types
% attention: dm consist of all particles with C++ types 1,2,3
PT_GAS = 1; PT_DM = 2; PT_STARS = 3; PT_BNDRY = 4;
HT_ALL = 1; HT_DM = 2; HT_GAS = 3; HT_STARS = 4; HT_BARY = 5;

sims = [ 'tarkin_13309_CSF_45x_snap_302';...
        'tarkin_21926_CSF_45x_snap_302';...
        'tarkin_25174_CSF_45x_snap_302'];
        
% temporary small snapshot 
% to be used for debugging purposes, isnap=3
%sims = ['tarkin_25174_CSF_10x_snap_009'];

dir_amd = '/data/achtland1/psteger/amd/';

prepare_m_arrays;

nsnap = size(sims); nsnap=nsnap(1);
for kaunt=1:nsnap
    sim = sims(kaunt,:);
    
    dir_output = [dir_amd 'output/' sim '/'];
    file_mtree = [dir_output 'ahf_out_mtree_idx2'];
    file_snap = [dir_output 'snapshot.hdf5'];
    file_halos = [dir_output 'ahf_out_halos2'];
    file_hpropi = [dir_output 'Hprop.txt'];
    file_hprop = [dir_output 'HpropStripped.txt'];
    file_tensor = [dir_output 'tensor_1.hdf5'];
    % we have to change the last number [0..10] 
    % if we want smoothing scale ..._3 * 0.5 Mpc/h of tidal field
    nhalo = importdata(file_hpropi);
    nhalo = nhalo(1,3);
    prepare_arrays;
    
    ahf_input;
    ahf_prepare;

    hp_input;
    hp_prepare;

    tf_input;
    tf_prepare;

    exc_prepare;
    exc_concat;
    
    ahf_concat;
    hp_concat;
    tf_concat;
    
end

%separate_for_environment;

% what time did we need?
toc
