function [ x_star, f_vals ] = MFISTA( b,A,At,x_start,lambda,L,iters,TV_Iters,x_true)
%MFISTA Implements MFISTA algorithm according to Amir Beck's paper:
% Solves the following convex problem:
%       x_star = arg min_x ||Ax-b||_2^2 + lambda*TV(x)
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
%   x_true   - the true solution to the problem for comparison purposes

global DebugFlag PlotFlag;

% The step size
t_k = 1;

% Initializing the temp variable y
y_k = x_start;
x_k = x_start;

% Computing the first fucntion value
A_y_k = A(y_k);
F_k = norm(A_y_k-b,'fro').^2 + 2*lambda(1)*TV(y_k);

% The function vals array
f_vals = F_k;

% Precomputing the projection A'b
At_b = At(b);

% We can use on fixed lambda value, or change its value throughout iterations
if isscalar(lambda)
    lambda_vec = ones(1,iters)*lambda;
else
    lambda_vec = linspace(lambda(1),lambda(2),iters);
end

% Resetting the psnr_vals vector, that holds the PSNR values during iterations
psnr_vals = psnr(x_start,x_true);

% Starting the waitbar
% wb = MyWaitbar(0,'Running MFISTA-TV...');
cpb = ConsoleProgressBar();



% Starting main iteration of MFISTA
cpb.start(); cpb.setText('FISTA-TV Algorithm');

% Opening a figure for debug and analysis purposes
if DebugFlag && PlotFlag 
    figure('Name','MFISTA Status','Position',[50 50 1280 720]);
end

for k =1:iters
    
    lambda_k = lambda_vec(k);
 
    % Original TV prox with gradient descent version of the paper by Amir
    z_k = FGP( y_k - 2/L * (At(A_y_k) - At_b ) , 2*lambda_k/L, TV_Iters);

% Set bilateral filter parameters.
    w     = 5;       % bilateral filter half-width
    sigma = [.2 6]; % bilateral filter standard deviations

% Apply bilateral filter to each image.
% bflt_img1 = bfilter2(img1,w,sigma);

%     z_k = bfilter2(Proj_C(z_k) ,w,sigma );
%         z_k = FGP(z_k,2*lambda_k/L,TV_Iters);

    
    
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
    
    % Computing the projection for next iteration
    A_y_k = A(y_k);

    F_k = norm(A_y_k-b,'fro')^2 + 2*lambda_k*TV(y_k);
    f_vals = [f_vals,F_k]; %#ok<AGROW>

    % Computing the PSNR compared to the ground truth
    psnr_vals = [psnr_vals,psnr(x_k,x_true)]; %#ok<AGROW>
    
    % If the psnr is going down after at least 10 iterations we stop
    if (psnr_vals(end) < psnr_vals(end-1))&&(k>10)
        break;
    end
    
    % Plotting iteration analysis if DebugFlag is on
    if DebugFlag && PlotFlag
        subplot(131);imagesc(x_k); colorbar('southoutside');
        subplot(132);loglog(f_vals,'LineWidth',2); title('Function value');
        grid on;
        subplot(133);semilogx(psnr_vals,'LineWidth',2);
        title(['PSNR value ',num2str(psnr_vals(end))]);
        grid on;
        suptitle(['MFISTA, Iteration ',num2str(k)]);
        drawnow;
    end
    
%     MyWaitbar(k/iters,wb);
    cpb.setValue(k/iters);
end
% close(wb);
cpb.stop();
x_star = x_k;
