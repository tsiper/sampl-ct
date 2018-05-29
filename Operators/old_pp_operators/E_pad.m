function [ I ] = E_pad( I,m)
%E Pads an image of size nxn to size mxn

[n1, n2] = size(I);
if (n1 ~= n2)
    error('E Operator - The image is not square');
end
n = n1;

d = m-n;
if mod(d,2)
    d = (d-1)/2;
    I = [zeros(n,d),I,zeros(n,d+1)];
else
    d = d/2;
    I = [zeros(n,d),I,zeros(n,d)];
end


end

