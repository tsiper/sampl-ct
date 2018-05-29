function [ p,q ] = Lt( x)
%LT Performs the derivative operation according to Amir Beck's paper.
% Basically returning the x and y derivatives for the respective p and q

p = -diff(x);
q = -diff(x.').';

end