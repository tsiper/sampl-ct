function [ n ] = tomo_radon_2D( A,lambda, dX, teta, nm, PGrid )

%teta is the rotation around the x axis

[~,N,L]= size(A);
A= lambda*A/(2*pi); %Phase->OPD

[y_i,x_i]= meshgrid((1:N) -(floor(N/2)+1), (1:N) -(floor(N/2)+1));
z_i= zeros(N);

F_SPACE= zeros(N,N,N);
avg_mat= zeros(N,N,N);
ones_vec= ones(N*N,1);

for p=1:L
rot_matrix_x= x_rotate3D(-teta(p));
loc_mat = round(rot_matrix_x*[x_i(:) y_i(:) z_i(:)]' + (floor(N/2)+1));
loc_mat(loc_mat>=N+1)=N;
loc_mat(loc_mat<=0)=1;
locInd = sub2ind([N N N], loc_mat(1,:)', loc_mat(2,:)', loc_mat(3,:)');

S= fftshift(fft2(ifftshift(A(:,:,p))))  *  (dX^2);

%% AT - Trying to trasnform only 1 line to see the radial pattern
S_oneline = fftshift(fft2(ifftshift(A(1,:,p))))  *  (dX^2);
figure;plot(PGrid(1:256,1),S_oneline); % PGrid(1:256,2),
%%

S_vec= S(:);
     
F_SPACE(locInd)= F_SPACE(locInd)+S_vec;
avg_mat(locInd)= avg_mat(locInd)+ones_vec;
end

avg_mat(avg_mat==0)=1; 
F_SPACE= F_SPACE./avg_mat; 

refractive_index= real((fftshift(ifftn(ifftshift(F_SPACE))) * ((1/dX)^3)));  

n = refractive_index+ nm*ones(size(refractive_index)); 