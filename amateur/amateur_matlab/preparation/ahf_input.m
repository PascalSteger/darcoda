% ahf_input: Get information from AHFstep

x=importdata(file_halos);

ahf_npart1=x(:,1); ahf_npart2=x(:,2);
ahf_pcx=x(:,3); ahf_pcy=x(:,4); ahf_pcz=x(:,5);
ahf_vcx=x(:,6); ahf_vcy=x(:,7); ahf_vcz=x(:,8);

ahf_mvir=x(:,9); ahf_rvir=x(:,10);
ahf_vmax=x(:,11);
ahf_rmax=x(:,12);
ahf_sigv=x(:,13);
ahf_lambda=x(:,14);

ahf_jx=x(:,15); ahf_jy=x(:,16); ahf_jz=x(:,17);

ahf_a=x(:,18); ahf_eax=x(:,19); ahf_eay=x(:,20); ahf_eaz=x(:,21);
ahf_b=x(:,22); ahf_ebx=x(:,23); ahf_eby=x(:,24); ahf_ebz=x(:,25);
ahf_c=x(:,26); ahf_ecx=x(:,27); ahf_ecy=x(:,28); ahf_ecz=x(:,29);
ahf_ovdens=x(:,30);
ahf_redge=x(:,31);
ahf_nbins=x(:,32);
ahf_ekin=x(:,33); ahf_epot=x(:,34);
ahf_mbp_offset=x(:,35);
ahf_com_offset=x(:,36);
ahf_r2=x(:,37);
ahf_lambdaE=x(:,38);
ahf_n_gas=x(:,39);
ahf_mvir_gas=x(:,40);
ahf_lambda_gas=x(:,41);
ahf_jx_gas=x(:,42); ahf_jy_gas=x(:,43); ahf_jz_gas=x(:,44);
ahf_a_gas=x(:,45); ahf_eax_gas=x(:,46); ahf_eay_gas=x(:,47); ahf_eaz_gas=x(:,48);
ahf_b_gas=x(:,49); ahf_ebx_gas=x(:,50); ahf_eby_gas=x(:,51); ahf_ebz_gas=x(:,52);
ahf_c_gas=x(:,53); ahf_ecx_gas=x(:,54); ahf_ecy_gas=x(:,55); ahf_ecz_gas=x(:,56);

ahf_ekin_gas=x(:,57); ahf_epot_gas=x(:,58);
ahf_lambdaE_gas=x(:,59);
ahf_phi0=x(:,60);

y=importdata(file_mtree);
ahf_subno=y(:,1); ahf_hostno=y(:,2);
