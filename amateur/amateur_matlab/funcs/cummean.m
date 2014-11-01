function L=cummean(y)
        % compute running mean of (N,1) vector
L = mean([[0 y];[y 0]]);
L = circshift(L,[0 -1]);
L = L(1:length(L)-2);