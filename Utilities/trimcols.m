function [ A_trim ] = trimcols( A , M)
%TRIPCOLS Trims the matrix A columns to a size M

m = size(A,1);

if (m<M)
    warning('The trimming size is bigger than the actual column length');
    A_trim = A;
    return;
end

if (m==M);
    warning('The trimming size is the same as the actual column size');
    A_trim = A;
    return;
end


% Calculating the pad size ammount
trim_size = m-M;

trim_top    = floor(trim_size/2);
trim_bottom = ceil(trim_size/2);

A_trim = A(trim_top+1:end-trim_bottom,:,:);

end

