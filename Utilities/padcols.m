function [ A_pad ] = padcols( A,M )
%PADCOLS Symmetrically pads the columns in A to size N

[m,n] = size(A);

if (m>M)
    warning('The padding size is smaller than the actual column length');
    A_pad = A;
    return;
end

if (m==M);
%     warning('The padding size is the same as the actual column size');
    A_pad = A;
    return;
end

% Calculating the pad size ammount
pad_size = M-m;

pad_top    = floor(pad_size/2);
pad_bottom = ceil(pad_size/2);

A_pad = [zeros(pad_top,n); A; zeros(pad_bottom,n)];

end

