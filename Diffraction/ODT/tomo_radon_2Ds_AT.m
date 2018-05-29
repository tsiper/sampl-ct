function [ n_spurs,n_fbp ] = tomo_radon_2Ds_AT( A,lambda, pix_CCD, Mag, teta, nm )

%%PARAMETERS
[~,N,L]= size(A);
dX=pix_CCD/Mag; %Transversial jumps [m]

A= lambda*A/(2*pi); %OPD
scaleFactor = 1;
% Building polar grid with the right cooaridabtes for the K space values
%N_rd = 2*ceil(norm(size(x0)-floor((size(x0)-1)/2)-1))+3;
PolarGrid = BuildPolarGrid(N,L)/scaleFactor;

ang = 270;
rot_matrix= [cosd(-ang) -1*sind(-ang) ; sind(-ang) cosd(-ang)]; %around z
PolarGridRot = (rot_matrix*PolarGrid.').';

% For each slice seperatly the 2D FT of the object is built
refractive_index_AT_SPURS = zeros(N,N,N);
refractive_index_AT_FBP = zeros(N,N,N);


for t=1:N % AT - was 1:N

for a=1:L %running on all angles per each slice

S(a,:)= fftshift(fft(ifftshift(A(t,:,a))))  *  (dX); %DC at center

%MAP ON LINE
% for k=1:N      
%     F_space(new_vec(1,k,a), new_vec(2,k,a))= F_space(new_vec(1,k,a), new_vec(2,k,a))+S(a,k);
%     avg_mat(new_vec(1,k,a), new_vec(2,k,a))= avg_mat(new_vec(1,k,a), new_vec(2,k,a))+1;
% end

end

% IN the end, after mapping all angles for the slice-  we average (functions as LPF)
%  avg_mat(avg_mat==0)=1; %preventing division by zero
%  F_space= F_space./avg_mat; 

%3rd stage is a 2D IFFT per each slice
St = S.';
S_vec = St(:);

%figure;imagesc(real(invF1(St))); % Sinogram of St

%% Gili original fft
%refractive_index(:,:,t)= (fftshift(ifft2(ifftshift(F_space))) * ((1/dX)^2));  

%% Running SPURS
[OutputImage,~]= Run_SPURS_Diffraction(S_vec,PolarGridRot);
OutputImage = OutputImage*((1/dX)^2);
OutputImage = OutputImage / max(OutputImage(:)); % Values normalization
OutputRI = fliplr(imrotate(OutputImage,270)); % Accomodating to Gili's output
refractive_index_AT_SPURS(:,:,t) = OutputRI;

%% Running iradon
Sinogram = squeeze(A(t,:,:));
x_fbp = iradon(Sinogram,teta,N);
x_fbp = x_fbp*((1/dX)^2);
x_fbp = x_fbp / max(x_fbp(:));
x_fbp = fliplr(imrotate(x_fbp,180));
refractive_index_AT_FBP(:,:,t) = x_fbp;

end

%refractive_index=real(refractive_index); %reparing numerical errors

% figure;
% imagesc(OutputRI);title('ifft with shift');
%adding the medium through which the sample beam went through in the
%sample- free interferogram

n_spurs=refractive_index_AT_SPURS+ nm*ones(size(refractive_index_AT_SPURS)); %make up for differnce.
n_fbp=refractive_index_AT_FBP+ nm*ones(size(refractive_index_AT_FBP)); %make up for differnce.

%since we do beam refrencing.

