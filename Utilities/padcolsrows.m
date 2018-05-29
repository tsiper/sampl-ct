function [ A_pad ] = padcolsrows( A, m,n )
%PADCOLSROWS Pads the matrix to size m x n

A_pad = padcols(padcols(A,m).',n).';

end

