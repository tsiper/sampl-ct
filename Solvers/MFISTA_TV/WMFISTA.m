function [ x_star, f_vals ] = wMFISTA( b,A,At,x0,lambda,L,iters,TV_Iters,Wx,Wy )
%MFISTA Implements MFISTA algorithm according to Amir Beck's paper:
% 

global DebugFlag PlotFlag;

% The step size
t_k = 1;

% Initializing the temp variable y
y_k = x0;
x_k = x0;

F = @(x) norm(A(x)-b,'fro').^2 + 2*lambda*TV(x);

f_vals = F(x0);

wb = MyWaitbar(0,'Running MFISTA-TV...');

% Opening a figure for debug and analysis purposes
if DebugFlag && PlotFlag 
    figure('Position',[50 50 1280 720]);
end

At_b = At(b);


% Starting main iteration of MFISTA
for k =1:iters
    % Original TV prox with gradient descent version of the paper by Amir
    z_k = wFGP( y_k - 2/L * (At(A(y_k)) - At_b ) , 2*lambda/L, TV_Iters,Wx,Wy);
    
    % Making sure that the algorithm is monotone by checking the function value
    F_z = F(z_k);
    x_k_prev = x_k;
    if F_z < f_vals(end)
        x_k = z_k;
        f_vals = [f_vals,F_z]; %#ok<AGROW>
    else
        f_vals = [f_vals,f_vals(end)]; %#ok<AGROW>
    end
    
    % Updating the step size
    t_k_next = (1+sqrt(1+4*t_k^2)) / 2;
    
    % The gradient descent step size for the y_k vector
    y_k = x_k + (t_k / t_k_next)*(z_k-x_k) + (t_k-1)/t_k_next * (x_k - x_k_prev);
    
    
    % Plotting iteration analysis if DebugFlag is on
    if DebugFlag && PlotFlag
        subplot(121);imagesc(x_k);
        subplot(122);loglog(f_vals,'LineWidth',2);
        grid on;
        suptitle(['wMFISTA, Iteration ',num2str(k)]);
    end
    
    MyWaitbar(k/iters,wb);
end
close(wb);
x_star = x_k;
