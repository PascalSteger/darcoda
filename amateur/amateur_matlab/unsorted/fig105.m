% aim: reproduce fig.2 of Pereira, Bryan and Gill 2008
start;

% major axis is given by vector from biggest eigenvalue
% ahf_a > ahf_b > ahf_c, so take ahf_eax, ahf_eay, ahf_eaz

angle_r_em = cosine( ahf_eax_gas, ahf_eay_gas, ahf_eaz_gas, ahf_dxc, ahf_dyc, ahf_dzc );

[n,x] = hist(abs(angle_r_em(exc)),10);

bar(x,n/sum(~isnan(angle_r_em(exc))),1);
xlabel('cos(r, I_{gas})'); ylabel('fraction');