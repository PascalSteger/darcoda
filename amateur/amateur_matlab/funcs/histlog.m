%histlog:
function r = histlog(M,nbin,varargin)
        exl = M>0;
        
        lx = log10(M(exl));
        [n,x] = hist(lx,nbin);
        y = log10(n/sum(exl));
        
        if nargin <= 2
                plotopts = 'r';
        else
                plotopts = varargin{1};
        end

        plot(x,y,plotopts);
        errorbar(x,...
                y,...
                log10(n)./sqrt(n),...
                [plotopts '.']);
        r = stat(exl);