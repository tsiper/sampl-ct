function [ x_star, f_vals ] = FISTA_TV( b,A,At,x_start,lambda,L,iters,TV_Iters,x_true)
%MFISTA Implements MFISTA algorithm according to Amir Beck's paper:
% Solves the following convex problem:
%       x_star = arg min_x 0.5*||Ax-b||_2^2 + lambda*TV(x)
%
% Input:
% ======
%   b  - is the measurement vector or matrix that we wish to solve
%   A  - A function handle that can compute A(x)
%   At - A function handle to compute the adjoint operator A*(x)
%   x_start - The initial guess vector/matrix for the solution of x
%   lambda  - The regularization constant, usually 1e-3 - 1e-4 works nicely
%   L       - The Lipschitz constant. Affects step size, which is 1/L
%   iters   - The number of iteration for the main loop (a couple of tens-hunderds should suffice)
%   TV_Iters - Number of iteration for the inner TV optimization problem loop (20-30 should work)
%   x_true   - (optional) the true solution to the problem for comparison purposes

global DebugFlag;

% The step size
t_k = 1;

% Initializing the temp variable y
y_k = x_start;
x_k = x_start;

% Computing the first fucntion value
A_y_k = A(y_k);
F_k = 1/2*norm(A_y_k-b,'fro').^2 + lambda*TV(y_k);

% Resetting the minimal value as the first value
f_min = F_k;

% The function vals array
f_vals = F_k;

% Precomputing the projection A'b
At_b = At(b);

% Resetting the psnr_vals vector, that holds the PSNR values during iterations
if exist('x_true','var')
    psnr_vals = psnr(x_start,x_true);
end

% Starting the waitbar
% wb = MyWaitbar(0,'Running MFISTA-TV...');
cpb = ConsoleProgressBar();

% Starting main iteration of MFISTA
cpb.start(); cpb.setText('FISTA-TV Algorithm');

% Opening a figure for debug and analysis purposes
if DebugFlag
    figure('Name','MFISTA Status','Position',[50 50 1280 720]);
end

best_psnr_iter = 0;
for k =1:iters
     
    % Original TV prox with gradient descent version of the paper by Amir
    z_k = FGP( y_k - (1/L) * (At(A_y_k) - At_b ) , lambda/L, TV_Iters);
    
    % Making sure that the algorithm is monotone by checking the function value
    % %     F_z = F_k(z_k);
    x_k_prev = x_k;
    % %    if F_z < f_vals(end)
    x_k = z_k;
    % %    else
    % %        f_vals = [f_vals,f_vals(end)]; %#ok<AGROW>
    % %    end
    
    % Updating the step size
    t_k_next = (1+sqrt(1+4*t_k^2)) / 2;
    
    % The gradient descent step size for the y_k vector
    y_k = x_k + (t_k / t_k_next)*(z_k-x_k) + (t_k-1)/t_k_next * (x_k - x_k_prev);
    t_k = t_k_next;
    % Computing the projection for next iteration
    A_y_k = A(y_k);

    F_k = 1/2*norm(A_y_k-b,'fro')^2 + lambda*TV(y_k);
    
    % Taking the best solution and saving it
    if (F_k < f_min)||(k==1)
        f_min = F_k;
        iter_min = k;
        if ~exist('x_true','var') 
            x_star = y_k;
        end
    end
    f_vals = [f_vals,F_k]; %#ok<AGROW>

    % Computing the PSNR compared to the ground truth
    if exist('x_true','var') 
        % Taking the best solution PSNR-wise
        this_iter_psnr  = psnr(y_k,x_true);
        if this_iter_psnr > max(psnr_vals)
            x_star = y_k;
            best_psnr_iter = k;
        end
        psnr_vals = [psnr_vals,this_iter_psnr]; %#ok<AGROW>
%         if (psnr_vals(end) < psnr_vals(end-1))&&(k>10)
%             break;
%         end
        
    end    
    % If the psnr is going down after at least 10 iterations we stop
    
    % Plotting iteration analysis if DebugFlag is on
    if DebugFlag
        subplot(131);imagesc(x_k); colorbar('southoutside');
        subplot(132);loglog(f_vals,'LineWidth',2); 
        title(['Function value, Best at Iter ',num2str(iter_min)]);
        grid on;
        if exist('x_true','var')
            subplot(133);semilogx(psnr_vals,'LineWidth',2);
            title(['PSNR value ',num2str(psnr_vals(end)),' Best at Iter ',num2str(best_psnr_iter)]);
            grid on;
        end       
        suptitle(['MFISTA, Iteration ',num2str(k)]);
        drawnow;
    end
    
%     MyWaitbar(k/iters,wb);
    cpb.setValue(k/iters);
end
% close(wb);
cpb.stop();
