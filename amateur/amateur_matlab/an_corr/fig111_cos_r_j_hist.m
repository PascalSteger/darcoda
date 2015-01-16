% equivalent to fig. 3 of van den Bosch et al. 2002
% or better fig. 10 Croft et al. 2008
clf
set(gca,'FontSize',15);

%hp_cpa = (m_cosrj(1,:,4)); %1,:,halo type

% no clear correlation visible

m_ex = m_exc_1;
hp_cpa = (m_cosjt(1,logical(m_ex),1)); %tabc,:,halo type
hp_cpa = (m_cosjt(1,:,1));
% no correlation either

%hp_cpa = abs(m_cosrt(1,:,1)); % 1,:,tabc % very strong correlation

%hp_cpa = abs(m_cosrea(1,:,3)); % 1,:,halo type

%hp_cpa = abs(m_cosjea(1,logical(m_exc_2),2)); % 1,:,halo type
%hp_cpa = hp_cpa(logical(m_exc_1));

nbin = 10;
xx = linspace(0,1,nbin+1);
[n,x] = hist(hp_cpa,nbin);

% for i=1:nbin
%     xmed(i) = median(hp_cpa(ibin==i));
%     xeu(i)  = prctile(hp_cpa(ibin==i),84);
%     xel(i)  = prctile(hp_cpa(ibin==i),16);
%     
%     ymed(i) = n(i);
%     yeu(i)  = n(i)+sqrt(n(i));
%     yel(i)  = n(i)-sqrt(n(i));
% end

y = n/sum(m_ex)*nbin;
ey = y./sqrt(n);

errorbar(x,y,ey);
%bar(x,n);
%bar(x,n,1);

xlabel('cos (r, J_{all})'); ylabel('fraction');