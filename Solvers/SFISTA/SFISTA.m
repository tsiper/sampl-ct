function [ x_est, F_stack ] = SFISTA( b, Aop, Aop_t, Dop, Dop_t, Params ,x_start,x0 )
%SFISTA - Smoothing based FISTA solver for the sparse analysis problem
%
% The Lagrangian we are trying to minimize
% x_est = arg min ||Aop*x-b||_2^2 + Params.Lambda*||Dop_t*x||_1
%              x
%
% Syntax:
% -------
% [ x_est ] = SFISTA( b, A, Dop, Lf, Params )
%
% Inputs:
% -------
% b      - M long input data vector
% Aop    - M X N sensing matrix, or operator
% Aop_t  - Adjoint of operator Aop. If Aop is a matrix use Aop_t = []
% Dop    - Decomposition matrix (usually a basis). Can be an operator
% Dop_t  - Adjoint of operator Dop. If Aop is a matrix use Dop_t = []
% Params - Algorithm parameters
%          Lambda - 
%          mu0 - Initial smoothing parameter
%          muf - Final smoothing parameter
%          gamma -
%          Continuation
%          L - The Lipschitz constant
%
%
% Outputs:
% --------
% x_est  - Estimated signal
%

%% Getting the globals
global DebugFlag PlotFlag;

%% Input checks
% -------------------------------------------------------------------------
% Sensing matrix 
if isa(Aop, 'function_handle')                  % A is an operator
    A  = Aop;
    At = Aop_t;
else                                            % A is a matrix
    A  = @(x) Aop*x;
    At = @(x) Aop'*x;
end

% Decomposition matrix
if isa(Dop, 'function_handle')                  % D is an operator
    D  = Dop;
    Dt = Dop_t;
else                                            % D is a matrix
    D  = @(x) Dop*x;
    Dt = @(x) Dop'*x;
end

% Optional continuation parameters
try mu0 = Params.mu0; catch 
    mu0 = 1e1; 
end
try muf = Params.muf; catch 
    muf = 1e-6; 
end
try gamma = Params.gamma; catch 
    gamma = 2; 
end    

%% Initialization
% -------------------------------------------------------------------------
% Continuation parameters
if Params.Continuation                                    % Do continuation
    mu    = mu0;
else                                                % Don't do continuation
    mu    = Params.mu;   
    muf   = 0.1*mu;
end
ContFlag  = 1;                                          % Continuation flag
Lambda    = Params.Lambda; 
Lambda_mu = Lambda*mu;
L         = Params.L;                                  % Lipschitz constant                  

% Function value initialization

% Initialize x - problem is convex, so we can initialize with every vector we want
if nargin < 7
    x = At(b);
else
    x = x_start;
end

F_stack = FuncValue(x, b, A, Dt(x), mu, Lambda);
snr_stack = [];

wb_sfista = MyWaitbar(0,'Running SFISTA...');


if DebugFlag && PlotFlag
    figure('Position',[50 50 1280 720]);
end


%% Continuation loop
% -------------------------------------------------------------------------
while mu > muf && ContFlag
    % Initialization for the continuation loop
    t = 1;
    y = x;
    
    % Interations - maximum number of iterations per each continuation loop
    % ---------------------------------------------------------------------
    
    for ii = 1:Params.IterMax
        % Testing adaptive lambda
        Lambda = Lambda*(1-1/Params.IterMax);
        
        % Calculate gradients
        Grad_F = At(A(y) - b);                                             % Gradient of f at the point x_{k-1}
        Dt_x   = Dt(x);
        Grad_G = D(Dt_x - Soft(Dt_x, Lambda*mu))/mu;                     % Gradient of g at the point D^*x_{k-1}
        
        % Calculate gradient step
        z = y - (Grad_F + Grad_G)/L;
        
        % "Smart" update
        t_prev = t;
        t = 0.5 + 0.5*sqrt(1 + 4*t^2);
        
        % Apply monotonicity (Monotone FISTA) - This part is not mandatory
        x_prev = x;                                                        % x_prev = x_{k-1}
        Fk = FuncValue(z, b, A, Dt_x, mu, Lambda);
        % If we improved in this iteration (monotonicity)
        if F_stack(end) > Fk     % x = argmin{H_mu(x): z_k, x_{k-1}}                                  
            x = z;
%             break;
            F_stack = [F_stack Fk]; %#ok<AGROW>
        else
            F_stack = [F_stack, F_stack(end)]; %#ok<AGROW>
        end
%         x = z;
        
        % Update function value stack
        
        
        % Update y
        y = x + (t_prev/t)*(z - x) + ((t_prev - 1)/t)*(x - x_prev);
        
        
        if DebugFlag && PlotFlag
            snr_stack = [snr_stack, psnr(x,x0)]; %#ok<AGROW>
            subplot(131);imagesc(x);
            subplot(132);loglog(F_stack,'LineWidth',2);
            title('Lagrangian Value')
            grid on;
            subplot(133);semilogx(snr_stack,'LineWidth',2);
            title('PSNR');            grid on;
            % PSNR plot
            suptitle(['SFISTA, Iteration ',num2str(ii)]);
        end
                
        MyWaitbar(ii/Params.IterMax,wb_sfista);
    end
    
    
    if Params.Continuation
        mu = mu/gamma;                                    % Update mu
       
    else
        ContFlag = 0;
    end
    
    display(mu);
    
    
end

    close(wb_sfista);


% Output 
x_est = x;

%% Auxiliary functions
% -------------------------------------------------------------------------
% Soft thresholding with parameter alpha
function y = Soft(z, alpha)
y = sign(z).*max(abs(z) - alpha, 0);

% Calculation of the smooth function H_{\mu}(x)
function f = FuncValue(x, b, A, Dt_x, mu, lambda)
f_val = 0.5*norm(A(x) - b, 2)^2;
g_val = lambda*sum(Huber(Dt_x(:), lambda*mu));

f = f_val + g_val;
