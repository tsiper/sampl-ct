%%
N = 32;
Mp = 32;
theta = (0:1/Mp:1-1/Mp)*180;

Nd = round(sqrt(2)*N);


A = paralleltomo(N,theta);

x0 = LoadPhantom(N,'brain');


Aop  = @(x) reshape(A*x(:),[Nd,Mp]); 
AopT = @(y) reshape(A'*y(:),[N N]);

y = Aop(x0);

%% Plotting
figure;
subplot(121);imshow(y,[]);
subplot(122);imagesc(x0);
% axis square;
colormap('bone');