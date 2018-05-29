function [ I ] = vec2im( v )
%VEC2IM Receives a vector of size N^2 and returns the image NxN

% Making sure it is a vector
v = v(:);
% Measuring the size
N = length(v);
n = ceil(sqrt(N));

if n ~= sqrt(N)
    warning('The vector does not depict a square image of size NxN');
    v = [v; zeros(n^2-N,1)];
end

I = reshape(v,n,n);
end

