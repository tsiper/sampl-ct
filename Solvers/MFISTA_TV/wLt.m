function [ p,q ] = wLt( x, Wx, Wy )
%LT Performs the derivative operation according to Amir Beck's paper.
% Basically returning the x and y derivatives for the respective p and q

p = -Wy*diff(x)/(Wx+Wy);
q = -Wx*diff(x.').' /(Wx+Wy);

end