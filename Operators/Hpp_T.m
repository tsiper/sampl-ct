function [ x ] = Hpp_T( y )
%HPP Creates a good operator for doing Hpp or Hpp_T so it complies with HtH
% An even sided square image is required
m = size(y,1);

% N = (m-1)/2;
N = (m-1)/4;

% Padding y to bigger size
% y     = padcols(y,4*N+1);
x_pad = Rpp_T(M(y));

% Trimming down to x
x = x_pad(N/2+1:3*N/2,N/2+1:3*N/2);


end

