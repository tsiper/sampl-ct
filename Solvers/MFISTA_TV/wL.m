function [ x ] = wL( p,q, Wx, Wy )
%L Summary of this function goes here
%   Detailed explanation goes here

% The dimensions of the objective
m = size(q,1);
n = size(p,2);

% First we pad with zeros so the dimensions are good
p = [zeros(1,n);p]; % Padding to (m+1)x(n)
q = [zeros(m,1),q]; % Padding to (m)x(n+1)

% Generating the matrices
p_tag = diff(p);
q_tag = diff(q,1,2);

% Padding again
p_tag = [p_tag;zeros(1,n)]; % Padding to (m)x(n)
q_tag = [q_tag,zeros(m,1)]; % Padding to mxn

% Initializing x
x = (p_tag*Wy + q_tag*Wx)/(Wy+Wx);

end

