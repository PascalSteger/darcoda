function figprep(varargin)
% prepare display with title, xlabel, ylabel

clf;
set(gca,'FontSize',15);
set(gca,'FontAngle','italic')
hold on;
grid on;
box on;
if nargin>0
        title(varargin{1});
        if nargin>1
                xlabel(varargin{2});
                ylabel(varargin{3});
        end
end