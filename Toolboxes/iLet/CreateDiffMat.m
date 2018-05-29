function [ Dh, Dv, D ] = CreateDiffMat( M, varargin )
%CREATEDIFFMAT Creates the horizontal and vertical differentiation matrices. Assume M X M images
% 
% Syntax:
% -------
% [ Dh, Dv, D ] = CreateDiffMat( M, BC_Type )
%
% Inputs:
% -------
% M       - Size of M X M matrix x
% BC_Type - Boundary condition type. Currently supports 'circular'
%
% Outputs:
% --------
% Dh - Horizontal M^2 X M^2 differentiation matrix
% Dv - Vertical M^2 X M^2 differentiation matrix
% D  - Concatenation of D = [Dh^T Dv^T]^T 
%
% Written by Oren Solomon.
%

%% Handle inputs
if nargin > 1
    BC_Type = varargin{1};
else
    BC_Type = 'circular';
end

%% Type of BC
switch lower(BC_Type)
    case 'reflexive'
        disp('Not supported yet.');
    case 'circular'
        % Horizontal
        tmp_Dh = zeros(M - 1, M^2);
        jj = 1;
        for ii = 1:M - 1
            tmp_Dh(ii, jj) = 1;
            tmp_Dh(ii, jj + M) = -1;
            
            jj = jj + M;
        end
        tmp_Dh(M, 1) = -1;
        tmp_Dh(M, M^2 - M + 1) = 1;
        
        kk = 0;
        Dh = [];
        for ii = 1:M:M^2
            Dh = [Dh; sparse(circshift(tmp_Dh, kk, 2))];
            
            kk = kk + 1;
        end
        
        % Alternative realization
%         kk = 0;
%         Dh = sparse(M^2);
%         for ii = 1:M:M^2
%             Dh(ii:ii + M - 1, :) = sparse(circshift(tmp_Dh, kk, 2));
%             kk = kk + 1;
%         end
        
        % -----------------------------------------------------------------
        % Vertical
        tmp_v       = eye(M) + diag(-ones(M - 1,1), 1);
        tmp_v(M, 1) = -1;
        
        Dv = sparse(M^2);
        for ii = 1:M:M^2
            Dv(ii:ii + M - 1, ii:ii + M - 1) = sparse(tmp_v);
        end
    otherwise
        error('CreateDiffMat: Unknown boundary conditions');
end

D = [Dh; Dv];