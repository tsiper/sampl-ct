function [ y ] = Hpp( x )
%HPP Creates a good operator for doing Hpp or Hpp_T so it complies with HtH
% An even sided square image is required
N = size(x,1);

% Padding x to bigger size
x_pad = padarray(x,[N/2,N/2]);

% Applying a bigger sized pseduo-polar transform
y = Rpp(x_pad);

% % Trimming the columns
% y = trimcols(y,2*N+1);


end

