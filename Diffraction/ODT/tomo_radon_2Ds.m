function [ n ] = tomo_radon_2Ds( A,lambda, pix_CCD, Mag, teta, nm )

%%PARAMETERS
[~,N,L]= size(A);
dX=pix_CCD/Mag; %Transversial jumps [m]

A= lambda*A/(2*pi); %OPD

%First, for all slices and angles we should build vectors:
ind_col= ((1:N) -(floor(N/2)+1))' ;%DC at center. for ex. for N= 256 we get [-128:127].
ind_row= zeros(size(ind_col));
new_vec= zeros(2, N, L-1);%"new_vec" is the length of s, and each of its values contains the location
%of the corresponding pixel in S. 

%Rotation vector
for a=1: L %running on all angles per each slice  
  rot_matrix= [cosd(-teta(a)) -1*sind(-teta(a)) ; sind(-teta(a)) cosd(-teta(a))]; %around z
  new_vec(:,:,a)= round(rot_matrix*([ind_row ind_col]')+ (floor(N/2)+1) ) ;%may have a problem with 0 and 257
end
new_vec(new_vec>=N+1)=N;
new_vec(new_vec<=0)=1;

% For each slice seperatly the 2D FT of the object is built
refractive_index= zeros(N,N,N);
n_start = 1; % Default - 1
n_end = N; % Default - N
for t=n_start:n_end
F_space= zeros(N);    
avg_mat= zeros(N);
S= zeros(L,N);

for a=1: L %running on all angles per each slice

S(a,:)= fftshift(fft(ifftshift(A(t,:,a))))  *  (dX); %DC at center

%MAP ON LINE
for k=1:N      
    F_space(new_vec(1,k,a), new_vec(2,k,a))= F_space(new_vec(1,k,a), new_vec(2,k,a))+S(a,k);
    avg_mat(new_vec(1,k,a), new_vec(2,k,a))= avg_mat(new_vec(1,k,a), new_vec(2,k,a))+1;
end

end
 %IN the end, after mapping all angles for the slice-  we average (functions as LPF)
 avg_mat(avg_mat==0)=1; %preventing division by zero
 F_space= F_space./avg_mat; 

 %3rd stage is a 2D IFFT per each slice
refractive_index(:,:,t)= (fftshift(ifft2(ifftshift(F_space))) * ((1/dX)^2));  
end

refractive_index=real(refractive_index); %reparing numerical errors
%figure;
%imagesc(refractive_index(:,:,128));title('ifft with shift - Gili');
%adding the medium through which the sample beam went through in the
%sample- free interferogram
n=refractive_index+ nm*ones(size(refractive_index)); %make up for differnce.
%since we do beam refrencing.

