function [ y ] = Proj_C( x )
%PROJ_C Projecting our solution to the convex space we're working on, in this
%case it is the 0..1 box

% The box parameters
MaxVal = 1;
MinVal = 0;

% Projecting onto the box
y = x;
y(x>MaxVal) = MaxVal;
y(x<MinVal) = MinVal;


end

