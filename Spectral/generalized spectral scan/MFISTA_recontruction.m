function [ a_hat,f_vals,psnr_vals] = MFISTA_recontruction( p_input,h,h_grad,a_input,...
    lambda,L,iters,TV_iters,RadonMaterials,PhantomRes,L_coeffs,toTV,toStop,base,w_true)
%MFISTA Implements MFISTA algorithm according to Amir Beck's paper:
% Solves the following convex problem:
%       x_star = arg min_x ||Ax-b||_2^2 + lambda*TV(x)
%
% Input:
% ======
%   p_input  - is the measurement vector or matrix that we wish to solve
%   h  - A function handle that can compute A(x)
%   h_grad - A function handle to compute the gradient of h
%   a_input - The initial guess vector/matrix for the solution of x
%   lambda  - The regularization constant, usually 1e-3 - 1e-4 works nicely
%   L       - The Lipschitz constant. Affects step size, which is 1/L
%   iters   - The number of iteration for the main loop (a couple of tens-hunderds should suffice)
%   TV_Iters - Number of iteration for the inner TV optimization problem loop (20-30 should work)
%   RadonMaterials   - the true solution to the problem for comparison purposes
%   PhantomRes   - Phantom resolution for pseudo polar reconstruction
%   L_coeffs - p metric scaling coefficients
%   toTV   - flag stating whether to perform TV or not

%% parameters
global DebugFlag PlotFlag;

Y = PhantomRes; %phantom resolution
I = length(RadonMaterials); %Number of materials
K = length(p_input); %Number of spectras
MyPP = @(x) real(invF1(App(x))); %Pseudo Polar transformation
MyAdjPP = @(y) real(App_T(M(F1(y))));
slop_iter = 500;

AppAdjApp = @(x) MyAdjPP(MyPP(x));
MyInversePP = @(y) conjgrad( AppAdjApp , MyAdjPP(y), zeros(Y,Y) ,6);

DataSize = size(a_input{1});
D = DecOperator(2,'uniform');
[ AtA ] = BuildAtA( Y ,2,'uniform' );
%% initialization
% The step size
t_k = 1;

% Initializing the temp variables
y_k = a_input;
x_k = a_input;

% Computing the first fucntion value
h_y_k = h(y_k);
TV_w = 0;
for m = 1:I
%     w = MyInversePP(a_input{m});
%     w_true = MyInversePP(RadonMaterials{m});
%     TV_w = TV_w + TV(w);
TV_w = 0;
end

% optional to calculate the cost function value from just the fidelity term
% or to add the regularization term.
F_k = norm([h_y_k{:}]-[p_input{:}])^2; %+2*lambda(1)*TV_w;

% The function vals array
f_vals = F_k;
min_f_val = f_vals;

% We can use on fixed lambda value, or change its value throughout iterations
if isscalar(lambda)
    lambda_vec = ones(1,iters)*lambda;
else
    lambda_vec = linspace(lambda(1),lambda(2),iters);
end

% Resetting the psnr_vals vector, that holds the PSNR values during iterations
psnr_vals =0; % psnr(w,w_true);

% Starting the waitbar
cpb = ConsoleProgressBar();

% Starting main iteration of MFISTA
cpb.start(); cpb.setText('FISTA-TV Algorithm');
% Opening a figure for debug and analysis purposes
if DebugFlag && PlotFlag
    f2 = figure('Name','MFISTA Status','Position',[50 50 1280 720]);
end

%% Iterations
% Starting main iteration of MFISTA
delta_p = zeros(K,DataSize(1),DataSize(2));
arg_k = cell(I,1);
z_k = cell(I,1);
radonSize = DataSize(1)*DataSize(2);
update_y =zeros(I,radonSize);
[thetaInd,tInd] = ind2sub(DataSize,1:radonSize);
w_k = cell(I,1);
for k =1:iters
    lambda_k = lambda_vec(k);
    h_grad_y = h_grad(y_k); % cell(I,K)
    
    %% Gradient step
    %calculate the distance of the solution from measurements
    for n = 1:K
        delta_p(n,:,:)=h_y_k{n}-p_input{n}; % size: K X DataSize
    end
    %goes over sinogram indexs and updates
    %     update_y =[];
    %     for jj = 1:DataSize(2)
    %         for ii = 1:DataSize(1)
    %             h_grad_mat = h_grad_y(:,:,ii,jj); % size: I X K
    %             update_y =[update_y , 2/L.*h_grad_mat *(L_coeffs'.*delta_p(:,ii,jj))];%in rowstack
    %         end
    %     end
    parfor ii = 1:radonSize
        h_grad_mat = h_grad_y(:,:,thetaInd(ii),tInd(ii)); % size: I X K
        update_y(:,ii) = 2/L.*h_grad_mat *(L_coeffs'.*delta_p(:,thetaInd(ii),tInd(ii)));%in rowstack
    end
    for n=1:I
        temp = y_k{n}-reshape(update_y(n,:),DataSize);
        % project to positive values
        temp(temp<0)=0;
        arg_k{n} = temp;
    end
    %% Original TV prox with gradient descent version of the paper by Amir\
    DebugFlag=0;
    if PlotFlag
        for m = 1:I
            %             % Reverse radon transform for the material m
            %             FGP_arg = MyInversePP(arg_k{m});
            %             % Finding the norm of the spatial image for material m
            %             range_k = range(vec(FGP_arg));
            %             min_k = min(vec(FGP_arg));
            %             arg_k_norm = (FGP_arg-min_k)/range_k;
            %             %perform FGP
            %             w_k = FGP(arg_k_norm, 2*lambda_k/L, TV_iters);
            %
            %             % Normalizing back the image to its original intensity
            %             w_k_norm = w_k*range_k+min_k;
            %
            %             % Going back to the sinogram material space (and projecting to non-negative)
            %             temp = MyPP(w_k_norm);
            %             temp(temp<0) = 0;
            %             z_k{m}= temp;
            b = real(App_T(D(M(F1(arg_k{m})))));
            % Trimming unnecessary zeros
            b = b(Y/2+1:3*Y/2,Y/2+1:3*Y/2);
            
            % Finding the optimal Lipshitz constant of our AtA system
            L_fista = PowerMethod(AtA,Y);
            
            %% Running Fista for solving
            x_tv = FISTA_TV_Spectral(b,AtA,zeros(Y),L_fista,TV_iters);
            w_k{m} = double(x_tv);
            x_pad = padarray(w_k{m},[Y/2,Y/2]);
            z_k{m}  = Rpp(x_pad);
            
        end
    end
    DebugFlag=1;
    %% Step update
    x_k_prev = x_k;
    if toTV
        x_k = z_k;
    else
        x_k = arg_k;
    end
    
    t_k_next = (1+sqrt(1+4*t_k^2))/2;
    for m=1:I
        y_k{m} = x_k{m} + (t_k-1)/t_k_next*(x_k{m}-x_k_prev{m});
    end
    t_k = t_k_next;
    %% f_vals & psnr
    h_y_k = h(y_k);
    TV_w = 0;
    for m = 1:I
%         w = MyInversePP(y_k{m}); %maybe needs to be normalized
        TV_w = TV_w + TV(w_k{m});
    end
    
    F_k = norm([h_y_k{:}]-[p_input{:}])^2; %+2*lambda(1)*TV_w;
    f_vals = [f_vals,F_k];
    
    %saves best solutions
    if f_vals(end)<=min_f_val
        min_f_val = f_vals(end);
        a_hat = x_k;
    end
    
    % Computing the PSNR compared to the ground truth
    psnr_vals = [psnr_vals,psnr(w_true{1},w_k{1})];
    %     psnr_vals = [psnr_vals,psnr(x_true,x_k)];
    
    % If the f_vals is rising after at least 10 iterations we stop
    if toStop && k>slop_iter+1
        if (f_vals(end-slop_iter) == min_f_val)
            break;
        end
    end
    %% plotting
    % Plotting iteration analysis if DebugFlag is on
    
    if DebugFlag && PlotFlag  && ~mod(k,3)
        figure(f2)
        subplot(231);imagesc(x_k{1}); colorbar('southoutside');
        subplot(232);imagesc(x_k{2}); colorbar('southoutside');
        %         subplot(233);imagesc(x_k{3}); colorbar('southoutside');
        subplot(234);loglog(f_vals,'LineWidth',2); title(['Function value',num2str(f_vals(end))]);
        grid on;
        subplot(235);semilogx(psnr_vals,'LineWidth',2);
        title(['PSNR value ',num2str(psnr_vals(end))]);
        grid on;
        
        
        %         alpha_hat(:,:,1) = MyInversePP(x_k{1})./base(1);
        %         alpha_hat(:,:,2) = MyInversePP(x_k{2})./base(2);
        
        %         alpha_hat(:,:,3) = zeros(size(alpha_hat(:,:,2)));% MyInversePP(x_k{3})./base(3);
        subplot(233);imagesc(w_k{1});
        %         alpha_hat = (alpha_hat-min(alpha_hat(:)))/range(alpha_hat(:));
        subplot(236);imagesc(w_k{2});
        %         colorbar('southoutside');
        suptitle(['MFISTA, Iteration ',num2str(k)]);
        drawnow;
    end
    cpb.setValue(k/iters);
end

cpb.stop();
