%start;
figure
hp_lx = hp_l(1,:,4)./sqrt(sum(hp_l(:,:,4).^2,1));
hp_ly = hp_l(2,:,4)./sqrt(sum(hp_l(:,:,4).^2,1));
hp_lz = hp_l(3,:,4)./sqrt(sum(hp_l(:,:,4).^2,1));

hp_lxe = hp_lx; ahf_lxe = ahf_lx;

plot(hp_lxe,ahf_lxe,'x')
hp_lambdae = log10(hp_lambda(:,1,1));
ahf_lambdae = log10(ahf_lambda(1,:));
plot(ahf_lambdae,hp_lambdae,'.')
%plot(((hp_ly(exc))),((ahf_ly(exc))),'.')
%plot(((hp_lz(exc))),((ahf_lz(exc))),'.')

exc3 = abs(hp_lx(exc)-ahf_lx(exc))>std(hp_lx(exc)-ahf_lx(exc));
hp_lxt = hp_lx(exc); ahf_lxt = ahf_lx(exc);
ahf_mvirt = ahf_mvir(exc);
ahf_rvirt = ahf_rvir(exc);
hist(ahf_rvirt(exc3))
axis square;