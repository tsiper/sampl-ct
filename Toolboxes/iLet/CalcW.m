function [ W ] = CalcW( M, x, varargin )
%CALCW calculates the weightimg matrix for the IRLS loop of TV_iLET.
% 
% Syntax:
% -------
% [ W ] = CalcW( M, x, BC_Type )
%
% Inputs:
% -------
% M       - Size of M X M matrix x
% x       - Input matrix (or vector) for the weights
% BC_Type - Boundary condition type. Currently supports 'circular'
%
% Output:
% -------
% W       - 2M^2 X 2M^2 weighting matrix
%
% Written by Oren Solomon.
%

%% Handle inputs
if nargin > 2
    BC_Type = varargin{1};
else
    BC_Type = 'circular';
end

if nargin > 3
    Dh = varargin{2};
    Dv = varargin{3};
else
    % Create difference matrices
    [ Dh, Dv ] = CreateDiffMat( M, BC_Type );
end
 
% Just in case
x = x(:);

%% Weights
dd = 1./sqrt((Dh*x).^2 + (Dv*x).^2 + 1e-10);

% Create weighting matrix
W = spdiags([dd; dd], 0, sparse(2*M^2, 2*M^2));


