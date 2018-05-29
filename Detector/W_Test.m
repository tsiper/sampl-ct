%% Playing around with some decimation matrices for LR and HR sinograms

%% Initializing
Initialize();
clear;
close all;
clc;
N = 64;
M = 40;



%% Loading the phantom
x0 = LoadPhantom(N,'zubal');
theta_hr = (0:1/M:1-1/M)*180;

A_radon = paralleltomo(N,theta_hr,N,N);

y0 = reshape(A_radon*x0(:),[N,M]);

%% Designing a blurring matrix
D = 0.5*diag(ones(N,1))+diag(ones(N-1,1),1)+diag(ones(N-2,1),2)+diag(ones(N-3,1),3);
D = D+D';
figure;
subplot(221);imagesc(x0);
subplot(222);imagesc(D*x0);
subplot(223);imagesc(x0*D);
subplot(224);imagesc(D*x0*D);

% The kronecker product for working with vectors
Dbig = kron(D,D);
Dbig = kron(eye(N),D);
x_blur = Dbig*x0(:);

figure;
imagesc(vec2im(x_blur));



%% Creating some matrices that delete values from the sinogram
W1 = diag(repmat([1 0],[1 (N*M)/2]));
W2 = diag(repmat([ones(1,N),zeros(1,N)],[1,M/2]));
W3 = diag(repmat([1 0 0],[1 ceil((N*M)/3)]));
while length(W3)>N*M
    W3(:,length(W3))=[];
    W3(length(W3),:)=[];
end
W5 = diag(repmat([0 0 1],[1 ceil((N*M)/3)]));
while length(W5)>N*M
    W5(:,length(W5))=[];
    W5(length(W5),:)=[];
end
W6 = diag(repmat([0 1],[1 ceil((N*M)/2)]));
while length(W6)>N*M
    W6(:,length(W6))=[];
    W6(length(W6),:)=[];
end
W7 = diag(repmat([zeros(1,N),ones(1,N)],[1,M/2]));
while length(W7)>N*M
    W7(:,length(W7))=[];
    W7(length(W7),:)=[];
end
W8=zeros(M*N,M*N);
for i=2:M*N
    for j=2:M*N
        W8(i,j)=W3(i-1,j-1);
    end
end

%% Deleting values from the sinograms
y_dec1  = reshape(W1*y0(:),[N,M]);
y_dec2 = reshape(W2*y0(:),[N,M]);
y_dec3 = reshape(W3*y0(:),[N,M]);
y_dec5 = reshape(W5*y0(:),[N,M]);
y_dec6 = reshape(W6*y0(:),[N,M]);
y_dec7 = reshape(W7*y0(:),[N,M]);
y_dec8 = reshape(W8*y0(:),[N,M]);

% add definitions for y_dec3,y_dec4 and so on

%% Plotting
figure(1);
subplot(331);imagesc(x0);
subplot(332);imagesc(y0);
subplot(333);imagesc(y_dec1);
subplot(334);imagesc(y_dec6);
subplot(335);imagesc(y_dec2);
subplot(336);imagesc(y_dec7);
subplot(337);imagesc(y_dec3);
subplot(338);imagesc(y_dec5);
subplot(339);imagesc(y_dec8);
% add subplots for y_dec3,y_dec4 and so on

%%
figure(2)
subplot(331);spy(W1);
subplot(332);spy(W6);
subplot(333);spy(W2);
subplot(334);spy(W7);
subplot(335);spy(W3);
subplot(336);spy(W5);
subplot(337);spy(W8);

% add spy for W3,W4 and so on


%% calculate HR

num_of_W=7; % change according the number of W matrices
size_w_m=length(y0(:));
size_w_n=size_w_m;
all_W=zeros(size_w_m*num_of_W,size_w_n);

all_LR_col=zeros(size_w_m*num_of_W,1);
all_W(1:size_w_m,1:size_w_n)=W1;
all_W(size_w_m+1:size_w_m*2,1:size_w_n)=W2;
all_W(size_w_m*2+1:size_w_m*3,1:size_w_n)=W3;
all_W(size_w_m*3+1:size_w_m*4,1:size_w_n)=W5;
all_W(size_w_m*4+1:size_w_m*5,1:size_w_n)=W6;
all_W(size_w_m*5+1:size_w_m*6,1:size_w_n)=W7;
all_W(size_w_m*6+1:size_w_m*7,1:size_w_n)=W8;


%% Just solving the big system


% add W3,W4 and so on
all_LR_col(1:size_w_m)=y_dec1(:);
all_LR_col(1+size_w_m:2*size_w_m)=y_dec2(:);
all_LR_col(1+2*size_w_m:3*size_w_m)=y_dec3(:);
all_LR_col(1+3*size_w_m:4*size_w_m)=y_dec5(:);
all_LR_col(1+4*size_w_m:5*size_w_m)=y_dec6(:);
all_LR_col(1+5*size_w_m:6*size_w_m)=y_dec7(:);
all_LR_col(1+6*size_w_m:7*size_w_m)=y_dec8(:);


%% Just plain solving
y_hr = all_W\all_LR_col;
y_hr = reshape(y_hr,[N M]);
% add y_dec3,y_dec4 and so on

%% Running fista
x_start=y_dec1(:);
field1 = 'L0';  value1 = 10*max(eig(all_W'*all_W));
field2 = 'eta';  value2 = 1.1;
field3 = 'Lambda';  value3 = 0;
field4 = 'IterMax';  value4 = 2000;
field5 = 'Type';  value5 = "const";
field6 = 'NonNegOrth';  value6 = 1;
Params = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6);
[ x_out, x_history ]=FISTA(all_W,all_LR_col,zeros(1280,1280),x_start(:),Params);
x_out_recon=reshape(x_out,N,M);

%% imshow the reconstruction sinogram and psnr calculate

figure(3);
subplot(121);imagesc(y0);
title('The Original Sinogram');
subplot(122);imagesc(x_out_recon);
title('The Reconstrected Sinogram');
Psnr=psnr(x_out_recon,y0);
