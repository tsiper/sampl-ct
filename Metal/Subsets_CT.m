x = LoadPhantom(256,'zubal');
y = GaussianNoise(radon(x),0.05);

y1 = y(:,1:3:end);
y2 = y(:,2:3:end);
y3 = y(:,3:3:end);


[U,S,V] = svd([y1(:)';y2(:)';y3(:)']);
% S(1,1) = 0;
% S(2,2) = 0;
S(3,3) = 0;

% Eliminating last rank
Y = U*(S*V);
y1 = Y(1,:)'; y1 = reshape(y1,[367,60]);
y2 = Y(2,:)'; y2 = reshape(y2,[367,60]);
y3 = Y(3,:)'; y3 = reshape(y3,[367,60]);

%%

x0 = iradon(y,0:179);
x1 = iradon(y1,0:3:179);
x2 = iradon(y2,1:3:179);
x3 = iradon(y3,2:3:179);





figure;
subplot(231);imagesc(x1);
subplot(232);imagesc(x2);
subplot(233);imagesc(x3);
subplot(234);imagesc(x1+x2+x3);
subplot(235);imagesc(x0);

x_sum = (x1+x2+x3)/3;
x_sum = x_sum(2:end-1,2:end-1);
x0    = x0(2:end-1,2:end-1);

%%
CompareImages(x_sum,x,'X-Sum');
CompareImages(x0,x,'X-Original');