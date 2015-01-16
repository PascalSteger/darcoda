function nn = nannan( someList )
% return someList stripped off the NaN values

nn = someList(not(isnan(someList)));
