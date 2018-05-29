function [ x, obj_stack, a_stack ] = TV_iLET( y, H, Params )
%TV_ILET performs TV iLET based reconstruction for the TV based problem
%
%    min_{x} 0.5||y - Hx||_2^2 + \lambda TV(x)
%
% This is not an efficient implementation!
% Syntax:
% -------
% [ x, obj_stack, a_stack ] = TV_iLET( y, H, Params )
% 
% Inputs:
% -------
% y      - Blurred image such that y = Hx
% H      - 2D convolution matrix of appropriate dimensions 
% Params - Algorithm parameters
%          BC          - Boundary conditions. Currently support 'circular' only
%          lambda      - Regularization parameter
%          SizeM       - Size M of the reconstructed M X M image x
%          L           - Lipschitz constant of the quadratic term
%          IterMax     - Number of iterations
%          IterMaxIRLS - Number of IRLS iterations
%
% Outputs:
% --------
% x         - SizeM X SizeM solution
% obj_stack - vector of the objective function evolution
% a_stack   - Matrix of the iLET coefficients evolution (each column corresponds to a different iteration)
%
% Writen by Oren Solomon.
%

global showValue updateReg VERBOSE

if VERBOSE; disp('Starting TV iLET deconvolution.'); end;
if VERBOSE; disp('-------------------------------'); end;

%% Parse inputs
BC_Type = Params.BC;
lambda  = Params.lambda;
M       = Params.SizeM; 

%% Initialization
% [M, N] = size(y);       % Assume M == N, x is the same size of y
y = y(:);

obj_stack = [];
a_stack   = [];

% Create convolution matrix from the PSF
% if VERBOSE; fprintf('Calculating H...'); end;
% H = convmtx2(PSF, M, N);
% if VERBOSE; disp('Done,'); end;

% Calculate the constants
if VERBOSE; fprintf('Calculating constants...'); end;
b    = H'*y;
% bF   = F'*b;
HTH  = H'*H;
% HTHF = F'*HTH*F;

if ~isempty(Params.L)
    L = Params.L;
else
    L = max(eig(HTH));
end

% Create difference matrices
[ Dh, Dv, D ] = CreateDiffMat( M, BC_Type );
% DF = D*F;
if VERBOSE; disp('Done.'); end

% Solution vector
mu     = [1 10]*L/1000;                  % Step size is 1/L
x      = (HTH + mu(1)*eye(M^2))\b;
x_prev = (HTH + mu(2)*eye(M^2))\b;

% Initialize F - check this again
idim = 5;
F    = zeros(M^2, idim);

% Initialize a
a    = ones(idim, 1)*1e-1; 

%% Iterations
for ii = 1:Params.IterMax
    % Build the LET basis
    % ---------------------------------------------------------------------
    grad_step = x - (1/L)*(HTH*x - b);  % All in vector form (step size is 1/L)
    
    % TV minimization base vector.
    
    % Prox calculation part
    % ---------------------------------------------
    % Denoising, based on Amir Beck's paper: "Fast Gradient-Based Algorithms for Constrained
    % Total  Variation  Image  Denoising and Deblurring Problems", IEEE TRANSACTIONS ON IMAGE PROCESSING, VOL. 18, NO. 11, NOVEMBER 2009. 
    % This code snippet is taken from the function "deblur_tv_fista.m". TV denoising includes box constraints such as non-negativity.
    
    % Projection onto the non-negative orthant, only for a non-negative constraint
    NonNegOrth = 0;
    if NonNegOrth == 1
        lb = 0;
    else
        lb = -Inf;
    end
    
    % Invoking the denoising procedure. Upper bound is Inf.
    parsin.MAXITER = 20;
    parsin.epsilon = 1e-5;
    parsin.print   = 0;
    parsin.tv      = 'iso';             % 'iso': Isotropic TV norm, 'l1': Anisotropic TV norm
    parsin.fig     = 0;
%     L              = 18*lambda^2;       % Lipschitz constant - 16*lambda^2 according to Amir's paper. But this has higher PSNR
    
    if ii == 1
%         [x_TV_tmp,~,~,P] = denoise_bound_init(reshape(grad_step, M, M), lambda/L, lb, Inf, [] ,parsin);           % First iteration
        [x_TV_tmp,~,~,P] = denoise_bound_init(reshape(grad_step, M, M), lambda/L, 0, 1, [], parsin);            % Other iterations with warm start
    else
%         [x_TV_tmp,~,~,P] = denoise_bound_init(reshape(grad_step, M, M), lambda/L, lb, Inf, P, parsin);            % Other iterations with warm start
        [x_TV_tmp,~,~,P] = denoise_bound_init(reshape(grad_step, M, M), lambda/L, 0, 1, P, parsin);            % Other iterations with warm start
    end
    
    % Update the basis matrix F
    % ---------------------------------------------------------------------
    F(:, 1) = x_prev;
    F(:, 2) = x;
    F(:, 3) = x_TV_tmp(:);                              % It is posible to take also x - x_TV_tmp(:)
    F(:, 4) = (HTH + mu(1)*eye(M^2))\x_TV_tmp(:);       % It is posible to take also (HTH + mu(1)*eye(M^2))\(x - x_TV_tmp(:))
    F(:, 5) = (HTH + mu(2)*eye(M^2))\x_TV_tmp(:);       % It is posible to take also (HTH + mu(2)*eye(M^2))\(x - x_TV_tmp(:))
%     F(:, 6) = (HTH + mu(3)*eye(M^2))\x_TV_tmp(:);       % It is posible to take also (HTH + mu(2)*eye(M^2))\(x - x_TV_tmp(:))
%     F(:, 7) = (HTH + mu(4)*eye(M^2))\x_TV_tmp(:);       % It is posible to take also (HTH + mu(2)*eye(M^2))\(x - x_TV_tmp(:))
%     
    % Update the LET weights by IRLS type minimization
    % ---------------------------------------------------------------------
    for jj = 1:Params.IterMaxIRLS
        W = CalcW( M, F*a , 'circular', Dh, Dv );
%         a = (HTHF + 2*lambda*DF'*W*DF)\bF;
%         a = (F'*(HTH + 2*lambda*D'*W*D)*F)\(F'*b);
            a = pinv(F'*(HTH + 2*lambda*D'*W*D)*F)*(F'*b);
    end

    % update the LET basis
    % ---------------------------------------------------------------------
    x_prev = x;
    x      = F*a; 
    
    a_stack = [a_stack a];
    
    % Store objective function values
    if showValue || updateReg
        obj = 0.5*norm(y - H*x, 2)^2 + lambda*TV_norm(Dh, Dv, x);
        obj_stack = [obj_stack obj];
    end   
    
    % Print to screen
    if showValue && VERBOSE
        fprintf('Iteration: %3d, Obj value: %15.5e\n', ii, obj_stack(end));
    end
    
end

% Output as an M X M image
x = reshape(x, M, M);

%% Auxiliary functions
function TVn = TV_norm(Dh, Dv, x)
TVn = sum(sqrt((Dh*x).^2 + (Dv*x).^2));




