function [ a_hat,f_vals,psnr_vals] = spectral_FISTA_reconstruction_OM( p_input,h,h_grad,a_input,...
    ReconstructionParams,w_true,invP_operator,P_operator)
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

Y = ReconstructionParams.PhantomRes; %phantom resolution
I = length(a_input); %Number of materials
K = length(p_input); %Number of spectras

DataSize = size(a_input{1});
D = DecOperator(2,'uniform');
[ HtH] = BuildAtA( Y ,2,'uniform' );
cg_max_iters = 100;
cg_err       = 1e-20;
%% initialization
% The step size
t_k = 1;

% Initializing the temp variables
y_k = a_input;
x_k = a_input;

% Computing the first fucntion value
h_y_k = h(y_k);
TV_w = 0;
w_true_toPSNR = [];
for m = 1:I
    w_true_toPSNR = [w_true_toPSNR; w_true{m}];
end

% optional to calculate the cost function value from just the fidelity term
% or to add the regularization term.
F_k = norm([h_y_k{:}]-[p_input{:}])^2;%+2*ReconstructionParams.lambda(1)*TV_w;

% The function vals array
f_vals = F_k;
min_f_val = f_vals;

% Resetting the psnr_vals vector, that holds the PSNR values during iterations
psnr_vals = 0;

% Starting the waitbar
cpb = ConsoleProgressBar();

% Starting main iteration of MFISTA
cpb.start(); cpb.setText('FISTA-TV Algorithm');

% Opening a figure for debug and analysis purposes
if DebugFlag && PlotFlag
    f2 = figure('Name','FISTA Status','Position',[50 50 1208 720]);% 800 500
end

%% Iterations
% Starting main iteration of FISTA
delta_p = zeros(K,DataSize(1),DataSize(2));
arg_k = cell(I,1);
z_k = cell(I,1);
radonSize = DataSize(1)*DataSize(2);
update_y =zeros(I,radonSize);
[thetaInd,tInd] = ind2sub(DataSize,1:radonSize);
w_k = cell(I,1);
a_mat = zeros(I,DataSize(1),DataSize(2));
temp_mat = zeros(I,radonSize);

if ReconstructionParams.toRecord 
    v = VideoWriter([ReconstructionParams.fileRecordName,'.avi']);
    open(v);
end
for k =1:ReconstructionParams.iters
    h_grad_y = h_grad(y_k); % cell(I,K)
    
    %% Gradient step
    %calculate the distance of the solution from measurements
    for n = 1:K
        delta_p(n,:,:)=h_y_k{n}-p_input{n}; % size: K X DataSize
    end
    %goes over sinogram indexs and updates
    parfor ii = 1:radonSize
        h_grad_mat = h_grad_y(:,:,thetaInd(ii),tInd(ii)); % size: I X K 
        update_y(:,ii) =2/ReconstructionParams.L.*diag(ReconstructionParams.La_coeffs)*(h_grad_mat *...
            (diag(ReconstructionParams.L_coeffs)*delta_p(:,thetaInd(ii),tInd(ii))));%in rowstack
    end
    norm_update = norm(update_y)
    if 1
        for n=1:I
            a_mat(n,:,:)=y_k{n}-reshape(update_y(n,:),DataSize);
        end
        % projection to positive values
        
        parfor ii = 1:radonSize
            temp_mat(:,ii)=ProjToPositive3_new(invP_operator,P_operator,a_mat(:,thetaInd(ii),tInd(ii)));
        end
        for n=1:I
            arg_k{n} = reshape(temp_mat(n,:),DataSize);
        end
    else
        for n=1:I
            arg_k{n} = y_k{n}-reshape(update_y(n,:),DataSize);
        end
    end
    %% Original TV prox with gradient descent version of the paper by Amir\
    w_toPSNR = [];
    DebugFlag = 0; %no plotting during FISTA
    for m = 1:I
        b = Hpp_T(D(arg_k{m}));
        temp_w = conjgrad(HtH,b,zeros(size(b)),cg_max_iters,cg_err);
        min_w = min(temp_w(:));
        range_w = range(temp_w(:));
        norm_w_k = (temp_w-min_w)/range_w;
        norm_x_k = FGP( norm_w_k , ReconstructionParams.lambda,...%*norm_update
            ReconstructionParams.TV_Iters);%
        w_k{m} = norm_x_k*range_w+min_w;
        w_toPSNR = [w_toPSNR; w_k{m}];
        z_k{m}  = Hpp(w_k{m});
    end
    DebugFlag = 1;
    %% Step update
    x_k_prev = x_k;
    if ReconstructionParams.toTV
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
        TV_w = TV_w + TV(w_k{m});
    end
    
    F_k = norm([h_y_k{:}]-[p_input{:}])^2;%+2*ReconstructionParams.lambda(1)*TV_w;
    f_vals = [f_vals,F_k];
    
    %saves best solutions
    if f_vals(end)<=min_f_val
        min_f_val = f_vals(end);
        a_hat = x_k;
    end
    
    % Computing the PSNR compared to the ground truth
    psnr_vals = [psnr_vals,psnr(w_true_toPSNR(:),w_toPSNR(:))];
    
    %% plotting
    % Plotting iteration analysis if DebugFlag is on
    
    if DebugFlag && PlotFlag  && ~mod(k,1)
        plotRecon3_OM(x_k,w_k,psnr_vals,f_vals,f2,k,invP_operator,ReconstructionParams.base);
        if ReconstructionParams.toRecord
            A = frame2im(getframe(f2));
            writeVideo(v,A);
        end
    end
    cpb.setValue(k/ReconstructionParams.iters);
end
if ReconstructionParams.toRecord 
    close(v);
end
cpb.stop();
