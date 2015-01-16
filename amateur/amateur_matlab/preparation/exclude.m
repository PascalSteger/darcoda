% exclude.m: definitions of excluded halos
% ahf and hp should already be prepared

% our starting point: all halos allowed
exc = ones(size(hp_mvir));

% hh is the relative amount of heavy DM particles
%hh = double(hp_ptypn(:,1,3)) + double(hp_ptypn(:,1,4));
%hh = hh./double(hp_ptypn(:,1,2));

% and should be smaller than a 1/100 (arbitrary),
%  the heavy particles should not contribute to more than this small amount
%exc = exc & hh' < 1e-2;              %sum(exc)

% and with a given minimal mass
% in units of 10^10/h MSun
exc = exc & hp_mvir > 0.01;  
%stat(exc)

% exclude halos too far away from cluster center
%exc = exc & dtot<10*ahf_rvir(1);     %sum(exc)

% only look at halos with more than a fixed number of particles (KNOB)
exc_1 = exc & hp_npart(:,1,1) >200;  % all
%stat(exc_1)
exc_2 = exc & hp_npart(:,1,2) >100;  % dm
exc_3 = exc & hp_npart(:,1,3) > 10;  % gas
exc_4 = exc & hp_npart(:,1,4) > 10;  % stars

% The above exclusions should be kept in mind when experiencing resolution effects.
