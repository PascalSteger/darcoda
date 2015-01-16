% ahf_prepare.m: process the raw data from ahf_input
% here only values that do not depend on hp_ are computed, for the latter, cf hp_prepare towards the end

ahf_cpa = cosine(ahf_jx, ahf_jy, ahf_jz, ahf_jx_gas, ahf_jy_gas, ahf_jz_gas);

ahf_dxc = ahf_pcx - ahf_pcx(1);
ahf_dyc = ahf_pcy - ahf_pcy(1);
ahf_dzc = ahf_pcz - ahf_pcz(1);

ahf_rtot = normN1(ahf_dxc, ahf_dyc, ahf_dzc);

