function [ rms_val ] = rmse( X,Y,mask)
%RMSE Summary of this function goes here
%   Detailed explanation goes here
if nargin==3
    X = X.*mask;
    Y = Y.*mask;
end
v = X(:)-Y(:);
n = length(v);
rms_val = sqrt(v'*v / n);

end

