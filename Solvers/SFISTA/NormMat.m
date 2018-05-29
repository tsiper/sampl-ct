function [ An, NrmMat ] = NormMat( A, NormType )
%NORMMAT - normalizes the columns of A, according to NormType
%
% Syntax:
% -------
% [ An, NrmMat ] = NormMat( A, NormType )
%
% Inputs:
% -------
% A        - Matrix to normalize
% NormType - Type of normalization to perform (currently supports only 'l2')
%
% Output:
% -------
% An       - Normalized (columns) matrix
% NrmMat   - Diagonal matrix of weights
%
% Ver 1: Written by Oren Solomon, Technion I.I.T, 04-10-2015
%

switch lower(NormType)
    case 'l2'
        NormACols = sqrt(diag(A'*A));
        NrmMat = diag(1./NormACols);
        An = A*NrmMat;
    otherwise
        error('NormMat: NormType not supported.');
end