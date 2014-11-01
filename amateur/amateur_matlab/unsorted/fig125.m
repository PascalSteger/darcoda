%start

nbins = 15;

%angle between stars and gas
hist(abs(m_cosphiauto(1,:,9)),nbins);

%angle between gas and DM
hist(abs(m_cosphiauto(1,:,10)),nbins);

%angle between stars and DM
%hist(abs(m_cosphiauto(1,:,16)),nbins);