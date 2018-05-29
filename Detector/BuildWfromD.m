function [ W ] = BuildWfromD( D,size_y )
%BUILDWFROMD D is the decc operator.
% size_y is the size of y [row,col] of the input sinogram of D
% W is the decimating matrix
M = size_y(1);
N = size_y(2);

W = sparse(N*M,N*M);
for i=1:N*M
    y = zeros(N*M,1);
    y(i) = 1;
    
    y_mtx  = reshape(y,[M,N]);
    W(:,i) = vec(D(y_mtx));
end


end

