function [  ] = SaveFigure( filename, varargin )
%SAVEFIGURE Saves the figure according to filename
global FIGURE_PATH SaveFlag;

% Default values for x and y of the figure
x = 1280;
y = 800;

% Checking inputs
if (nargin==2)
    disp('Please enter both X and Y coordinates');
end
if (nargin>2)
    x = varargin{1};
    y = varargin{2};
end
if (nargin>3)
    flag = varargin{3};
else
    flag = [];
end

if FIGURE_PATH
    h = gcf;
    set(h,'Name',filename,'Color',[1 1 1]);

    if SaveFlag
        % First we're saving the figure
        savefig([FIGURE_PATH,'/',filename]);
        % Now saving .png and .pdf files
        
        set(h,'Units','Inches');
        pos = get(h,'Position');
        set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        
        print('-dpng',[FIGURE_PATH,'/', filename]);
        print(h,[FIGURE_PATH,'/',filename],'-dpdf','-r0');
        close;
    end

end

end