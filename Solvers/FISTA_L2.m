function [ x_star ] = FISTA_L2( b, A, At, x_start, lambda, L, iters )
%FISTA_L2 Solves the minimization problem of:
% x_star = arg min_x ||A*x-y||_2^2 + lambda*||x||_2^2
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

global DebugFlag;
theta=floor((0:1/(177-1):1)*179);

% The step size
t_k = 1;

% Initializing the temp variable y
y_k = x_start;
x_k = x_start;

% Computing the first fucntion value
A_y_k = A(y_k);
F_k = norm(A_y_k-b,'fro')^2 + lambda*norm(y_k)^2;

% The function vals array
f_vals = F_k;
f_min  = F_k;

% Precomputing the projection A'b
At_b = At(b);

% Starting the waitbar
% cpb = ConsoleProgressBar();
% Starting main iteration of MFISTA
% cpb.start(); cpb.setText('FISTA-L2 Algorithm');

% Opening a figure for debug and analysis purposes
if DebugFlag
    figure('Name','MFISTA Status','Position',[50 50 1280 720]);
end

iter_min = 1;
for k =1:iters
     
    % Original TV prox with gradient descent version of the paper by Amir
    z_k = y_k - (2/L) * (At(A_y_k) - At_b );
    
    % The prox solution of the L2 norm
    z_k = z_k / (1+lambda);
%     z_k = wthresh(z_k,'s',lambda);

    x_k_prev = x_k;
    x_k = z_k;

    % Updating the step size
    t_k_next = (1+sqrt(1+4*t_k^2)) / 2;
    
    % The gradient descent step size for the y_k vector
    y_k = x_k + (t_k / t_k_next)*(z_k-x_k) + (t_k-1)/t_k_next * (x_k - x_k_prev);
    t_k = t_k_next;
    % Computing the projection for next iteration
    A_y_k = A(y_k);

    F_k = norm(A_y_k-b,'fro')^2 + lambda*norm(y_k)^2;
    
    % Taking the best solution and saving it
    if (F_k < f_min)||(k==1)
        f_min = F_k;
        iter_min = k;
        x_star = y_k;
    end
    
    f_vals = [f_vals,F_k]; %#ok<AGROW>
    
    % Plotting iteration analysis if DebugFlag is on
    if DebugFlag
        subplot(121);imshow(iradon(reshape(x_k,367,177),theta)); colorbar('southoutside');
        subplot(122);loglog(f_vals,'LineWidth',2); title(['Function value - min achieved at ',num2str(iter_min)]);
        grid on;
        suptitle(['FISTA-l2, Iteration ',num2str(k)]);
        drawnow;
    end
    
%     MyWaitbar(k/iters,wb);
%     cpb.setValue(k/iters);
%     cpb.setText(['FISTA-L2, Func Value= ',num2str(F_k)]);
end
% close(wb);
% cpb.stop();
