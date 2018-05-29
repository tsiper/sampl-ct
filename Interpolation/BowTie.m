function [ H ] = BowTie( N,M,delta )
%BOWTIE This function builds the bow-tie filter
% N - The number of detectors
% M - The number of projection angles

if nargin == 2 
    delta = 0.1;
end

% Calculating the sub-indices of the filter
m = linspace(-1,1,M);
n = linspace(-1,1,N);

[xx,yy] = meshgrid(m,n);

H = double((abs(yy)>=abs(xx)-delta));

% A little blurring
% avgH = integralKernel([1 1 7 7], 1/49);
% H = integralFilter(H, avgH);

end
