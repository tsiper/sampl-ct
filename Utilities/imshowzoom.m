function [ ] = imshowzoom( I ,Letter)
%IMSHOWZOOM Shows a zoom magnification of an image
    N = size(I,1);
%     I = imresize(I,[256 256]);

    TitleParams = {'FontSize',46,'FontWeight','bold','FontName','Times New Roman', ...
    'Color',[1 1 0],'Units','normalized'};

    %%  Paramters for the zoom region
    y = round(0.63*N);
    x = round(0.32*N);
%     y = 126;
%     x = 136;
    RegionSize = round(N/8);
    Factor = 3;
    Color = [.9 0 0];
    PosX = .025;
    PosY = .9;
    
    %% Designing and building the plot zoom
    
%     ZoomRegion = [x:x+RegionSize-1,y:y+RegionSize-1];
    ZoomImage = imresize(I(y:y+RegionSize-1,x:x+RegionSize-1),Factor,'nearest');
    G = I;
    
    % Contrast Enhancement for zoomed image
%     ZoomImage = imadjust(ZoomImage,[0.5 1],[0 1]);
   
%     DisplayRegion = [N-Factor*N+1:N,N-Factor*N+1:N];
    G(N-Factor*RegionSize+1:N,N-Factor*RegionSize+1:N) = ZoomImage;

    imagesc(G);axis off square;
%     colormap('default');
    colormap('bone');
    rectangle('Position',[x,y,RegionSize,RegionSize],'EdgeColor',Color,'LineWidth',3);
    rectangle('Position',...
        [N-RegionSize*Factor+1,N-RegionSize*Factor+1,RegionSize*Factor-1,RegionSize*Factor-1], ...
        'EdgeColor',Color,'LineWidth',3);
    line([x+RegionSize-1,N],[y,N-RegionSize*Factor],'LineWidth',3,'Color',Color);
    line([x,N-RegionSize*Factor],[y+RegionSize-1,N],'LineWidth',3,'Color',Color);

    if nargin == 2
        text(PosX,PosY,['(',Letter,')'],TitleParams{:});
    end

end

