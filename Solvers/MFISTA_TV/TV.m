function [ TV_Norm, D ] = TV( x )
%TV Summary of this function goes here
%   Detailed explanation goes here

[m,n] = size(x);

p = -diff(x);
q = -diff(x.').';


% Padding with zeros
p = [p;zeros(1,n)];
q = [q,zeros(m,1)];


% The TV matrix
D = sqrt(p.^2 + q.^2);

TV_Norm = sum(D(:));

end

