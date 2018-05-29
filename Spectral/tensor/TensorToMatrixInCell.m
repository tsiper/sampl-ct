function [ MatrixesInCells ] = TensorToMatrixInCell( Tensor , number_of_cells )
% Convert Tensor to matrixes in cells. 
    MatrixesInCells = cell(number_of_cells , 1);
    
    for i=1:size(Tensor,3)  %size is from the Tensor ToolBox. function m = size(tensor,idx)
        MatrixesInCells{i,1} = double(Tensor(:,:,i));
    end

end