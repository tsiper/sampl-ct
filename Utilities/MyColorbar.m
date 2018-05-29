function [  ] = MyColorbar( axis_handle ,varargin )
% Plots a colorbar to the right of the figure, without affecting its position
% and scale. axis_handle can also be a vector of axis handles. 

%You can specify additional colorbar options by using the varargin 

if nargin<1
    axis_handle = gca;
end

for i=1:length(axis_handle)
    axPos = axis_handle(i).Position;
    barPos = [axPos(1)+axPos(3)+0.02, axPos(2), 0.02, axPos(4)];
    colorbar(axis_handle(i),'location','manual','position',barPos,varargin{:});
end


end

