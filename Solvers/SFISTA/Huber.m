function [ y ] = Huber( x, alpha )
%HUBER - computes the output of the Huber function:%
%   y = (1/2/alpha)*(x^2) if |x| < alpha
%   y = |x| - alpha/2     else
%
% Syntax:
% -------
%       [ y ] = Huber( x, alpha )
%
% Inputs:
% -------
% x     - Input vector to compute the Huber function on.
% alpha - alpha > 0, a scalar parameter.
%
% Output:
% -------
% y     - Output vector with the same length as x
% 
% Ver. 1. Written by Oren Solomon, Technion I.I.T. 15-05-2016
%

% Assume x is a vector - always a column vector
x = x(:);

x_abs = abs(x);

y = (1/2/alpha)*(x.^2).*double(x_abs < alpha) + (x_abs - alpha/2).*double(x_abs >= alpha);


