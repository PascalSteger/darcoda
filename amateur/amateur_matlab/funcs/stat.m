function [mean,stdev,occupied] = stat(V)
        % finds mean, standard deviation and occupation fraction of a vector
        n = length(V);
        mean = sum(V)/n;
        stdev = sqrt(sum((V-mean).^2/n));
        occupied = sum(V~=0)/n;
