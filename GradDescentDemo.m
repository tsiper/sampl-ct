%% Initialize the demo
Initialize;

%% Create the CT System
N = 256;
sigma = 0.03;
f0 = LoadPhantom(N,'brain');
[A,~,~,theta,M,K] = paralleltomo(N);

%% Calculating
SparseMeasure = @(x) nnz(x)/numel(x);
SystemSparsity = SparseMeasure(A)*100;

%% The projections
p = A*f0(:);
p_n = p+sigma*randn(size(p))*range(p);

%% BackProjection
f_bp = A'*p_n;
f_bp = vec2im(f_bp);

%% Filtered Back Projection
f_fbp = fbp(A,p_n,theta);
f_fbp = vec2im(f_fbp');

%% Gradient Descent
f_k = zeros(size(f0(:)));
% Atp = A'*p_n;
mu = 2e-5;
figure;
iters = 50;
ErrPlot = [];
for k = 1:iters
    f_k = f_k - mu*A'*(A*f_k-p_n);
    subplot(121);imagesc(vec2im(f_k));colormap('bone');
    ErrPlot = [ErrPlot,norm(A*f_k-p,'fro')]; %#ok<AGROW>
    title(['Iter ',num2str(k)]);
    subplot(122);loglog(ErrPlot);grid on;
    drawnow;
    pause(0.01);
end
    
%% Plotting
figure;
subplot(231);imagesc(f0);title('Clean Phantom');
subplot(232);imagesc(f_fbp);title('Filtered Backprojection');
subplot(233);imagesc(f_bp);title('Backprojection');
subplot(234);imagesc(reshape(p,[M,length(theta)]));title('Sinogram');
subplot(235);imagesc(reshape(p_n,[M,length(theta)]));title('Noisy Sinogram');
subplot(236);imagesc(vec2im(f_k)); title('GD Solution');
colormap('bone');

