function [ x, obj_stack, a_stack ] = TV_iLET_CT( y, H, Ht, Params,x0 )
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

global showValue updateReg VERBOSE DebugFlag;

if VERBOSE; disp('Starting TV iLET deconvolution.'); end;
if VERBOSE; disp('-------------------------------'); end;

%% Parse inputs
lambda  = Params.lambda;
N       = size(x0,1);
L       = Params.L;
mu      = 1e-3;

%% Initialization
% [M, N] = size(y);       % Assume M == N, x is the same size of y
% y = y(:);

obj_stack = [];
a_stack   = [];

% invHTH = @(x) Mnew(vec2im(x))

% Create convolution matrix from the PSF
% if VERBOSE; fprintf('Calculating H...'); end;
% H = convmtx2(PSF, M, N);
% if VERBOSE; disp('Done,'); end;

% Calculate the constants
if VERBOSE; fprintf('Calculating constants...'); end;
b    = Ht(y);
% bF   = F'*b;
HtH  = @(x) Ht(H(x));
% HTHF = F'*HTH*F;
% 
% if ~isempty(Params.L)
%     L = Params.L;
% else
%     L = max(eig(HtH));
% end

% Create difference matrices
[ Dh, Dv, D ] = CreateDiffMat( N );
% DF = D*F;
if VERBOSE; disp('Done.'); end

% Solution vector
% mu     = [1 10]*L/1000;                  % Step size is 1/L
% x      = (HTH + mu(1)*eye(M^2))\b;
% x_prev = (HTH + mu(2)*eye(M^2))\b;
% x      = (HtH + mu(1)*eye(M^2))\b;
% x_prev = (HtH + mu(2)*eye(M^2))\b;
x = b/L;
x_prev = b/L;
% x = zeros(size(x0));
% x_prev = x;

% Initialize F - check this again
idim = 7;
F    = zeros(N^2, idim);

% Initialize a
a    = ones(idim, 1)*1e-1; 

if DebugFlag
    iLETfig = figure('Name','iLET Analysis','Position',[50 50 1280 720]);
    psnr_stack = [];
    l2_norm_stack = [];
    IRLS_Fig = figure();

end

% Starting the time counter
cpb = ConsoleProgressBar();
% Starting main iteration of MFISTA
cpb.start(); cpb.setText(['TV-iLET, lambda=',num2str(lambda)]);


%% Iterations
for ii = 1:Params.IterMax
    % Build the LET basis
    % ---------------------------------------------------------------------
    HtH_x = HtH(x);
    grad_step = x - (1/L)*(HtH_x - b);  % All in vector form (step size is 1/L)
    
    % Updating the Lambda value
    lambda = norm(HtH_x - b,'fro')/N/3;

    % Invoking the denoising procedure. Upper bound is Inf.
    parsin.MAXITER = 20;
    parsin.epsilon = 1e-5;
    parsin.print   = 0;
    parsin.tv      = 'iso';             % 'iso': Isotropic TV norm, 'l1': Anisotropic TV norm
    parsin.fig     = 0;
%     L              = 18*lambda^2;       % Lipschitz constant - 16*lambda^2 according to Amir's paper. But this has higher PSNR
    
%     if ii == 1
% %         [x_TV_tmp,~,~,P] = denoise_bound_init(reshape(grad_step, M, M), lambda/L, lb, Inf, [] ,parsin);           % First iteration
%         [x_TV_tmp,~,~,P] = denoise_bound_init(reshape(grad_step, M, M), lambda/L, 0, inf, [], parsin);            % Other iterations with warm start
%     else
% %         [x_TV_tmp,~,~,P] = denoise_bound_init(reshape(grad_step, M, M), lambda/L, lb, Inf, P, parsin);            % Other iterations with warm start
%         [x_TV_tmp,~,~,P] = denoise_bound_init(reshape(grad_step, M, M), lambda/L, 0, inf, P, parsin);            % Other iterations with warm start
%     end
    x_TV_tmp  = FGP(grad_step,lambda/L,parsin.MAXITER);
    x_TV_tmp2 = FGP(grad_step,lambda/L*3,parsin.MAXITER);
    x_TV_tmp3 = FGP(grad_step,lambda/L/3,parsin.MAXITER);
    x_TV_tmp4 = FGP(grad_step,lambda/L*6,parsin.MAXITER);
    x_TV_tmp5 = FGP(grad_step,lambda/L/6,parsin.MAXITER);
% %     
%     x_TV_tmp  = FGP_unconstrained(grad_step,lambda/L,parsin.MAXITER);
%     x_TV_tmp2 = FGP_unconstrained(grad_step,lambda/L*10,parsin.MAXITER);
%     x_TV_tmp3 = FGP_unconstrained(grad_step,lambda/L/10,parsin.MAXITER);
%     
%     x_TV_tmp  = FGP_unconstrained(grad_step,lambda,parsin.MAXITER);
%     x_TV_tmp2 = FGP_unconstrained(grad_step,lambda*10,parsin.MAXITER);
%     x_TV_tmp3 = FGP_unconstrained(grad_step,lambda/10,parsin.MAXITER);
%     
    % Update the basis matrix F
    % ---------------------------------------------------------------------
    F(:, 1) = x_prev(:);
    F(:, 2) = Proj_C(grad_step(:));
    F(:, 3) = x(:);

    F(:, 4) = x_TV_tmp(:);                              % It is posible to take also x - x_TV_tmp(:)
    F(:, 5) = x_TV_tmp2(:);
    F(:, 6) = x_TV_tmp3(:);
    F(:, 6) = x_TV_tmp4(:);
    F(:, 7) = x_TV_tmp5(:);

    % Update the LET weights by IRLS type minimization
    % ---------------------------------------------------------------------
    % Calculating the multiplication from the left
    HtH_F = zeros(size(F));
    
    for n = 1:size(F,2)
        HtH_F(:,n) = vec(HtH(vec2im(F(:,n))));
    end
    
    for jj = 1:Params.IterMaxIRLS
        W = CalcW( N, F*a , 'circular', Dh, Dv );
%         a = (HTHF + 2*lambda*DF'*W*DF)\bF;
%         a = (F'*(HTH + 2*lambda*D'*W*D)*F)\(F'*b);
        DWD_F = zeros(size(F));
        for n = 1:size(F,2)
            DWD_F(:,n) = D'*(W*(D*F(:,n)));
        end
%         a = (F'*(HtH_F+2*lambda*DWD_F)+1e-1*eye(size(F,2)))\(F'*vec(b));
%         a = pinv(F'*(HtH_F+2*lambda*DWD_F)+mu*eye(size(F,2)))*(F'*vec(b));
%         a = (F'*(HtH_F+2*lambda*DWD_F)+mu*eye(size(F,2)))\(F'*vec(b));
          A_cg = F'*(HtH_F+2*lambda*DWD_F);
          b_cg = F'*vec(b);
%           a = pcg(A_cg,b_cg);
          cg_iters = 10;
          a = conjgrad(@(b) A_cg*b,b_cg,a,cg_iters);
        
%         if DebugFlag 
%             figure(IRLS_Fig);
%             imagesc(vec2im(F*a));
%             drawnow;
%         end
%         a = (F'*(HtH_F+2*lambda*DWD_F))\(F'*vec(b));
%             a = pinv(F'*(HtH + 2*lambda*D'*W*D)*F)*(F'*b);
    end

    % update the LET basis
    % ---------------------------------------------------------------------
    x_prev = x;
    x      = Proj_C(vec2im(F*a)); 
    
    a_stack = [a_stack a]; %#ok<AGROW>
    
    % Store objective function values
    if showValue || updateReg
        obj = 0.5*norm(y - H(x), 'fro')^2 + lambda*TV_norm(Dh, Dv, x);
        obj_stack = [obj_stack obj]; %#ok<AGROW>
    end   
    
    % Print to screen
    if showValue && VERBOSE
        fprintf('Iteration: %3d, Obj value: %15.5e\n', ii, obj_stack(end));
    end
    
    if DebugFlag
        
        obj = 0.5*norm(y - H(x), 'fro')^2 + lambda*TV_norm(Dh, Dv, x);
        obj_stack = [obj_stack obj]; %#ok<AGROW>
        
        psnr_val  = psnr(x,x0);
        psnr_stack = [psnr_stack,psnr_val]; %#ok<AGROW>
        
        l2_norm = 0.5*norm(HtH_x-b, 'fro')^2;
        l2_norm_stack = [l2_norm_stack,l2_norm]; %#ok<AGROW>
        
        figure(iLETfig);
        subplot(131);imagesc(x);
        colorbar('southoutside');
        subplot(232);semilogy(obj_stack); grid on;
        title(['Objective Value=',num2str(obj)]);
        subplot(235);semilogy(l2_norm_stack); grid on;
        title(['L2 norm Value=',num2str(l2_norm)]);
        subplot(133);plot(psnr_stack); grid on;
        title(['PSNR, Val= ',num2str(psnr_val),'dB']);
        drawnow;
    end
    
    cpb.setText(['TV-iLET, lambda=',num2str(lambda)]);
    cpb.setValue(ii/Params.IterMax);

end
cpb.stop();

%% Auxiliary functions
function TVn = TV_norm(Dh, Dv, x)
TVn = sum(sqrt((Dh*x(:)).^2 + (Dv*x(:)).^2));




