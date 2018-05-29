Initialize;
N = 128;

x = LoadPhantom(N,'zubal');
theta = 0:2:179;

Amat = paralleltomo(N,theta,N,N);

A  = @(x) reshape(Amat*x(:),[N,length(theta)]);
At = @(y) reshape(Amat'*y(:),[N N]);

y = A(x);
x_hat = At(y);
x_fbp = reshape( fbp(Amat,y(:),theta), [N N]);

figure;
subplot(221);imagesc(x);
subplot(222);imagesc(y);
subplot(223);imagesc(x_hat);
subplot(224);imagesc(x_fbp);

