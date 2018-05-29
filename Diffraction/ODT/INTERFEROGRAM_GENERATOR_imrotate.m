function [ I ] = INTERFEROGRAM_GENERATOR_imrotate( n, teta, lambda, alpha, medium, dX, pix_ccd )

%Preparations
N= length(teta); %Number of interferograms. 
[size_x,size_y,~]= size(n);
I= zeros(size_x, size_y, N+1);%preparing the array of interferograms

%Reference beam
%Reference beam amplitude
R= ones(size_x, size_y); %Refrence beam. common phase such as kz-wt is ignored since in substraction the common phase
R_n= ones(size(n)); %goes through air instead of sample. 
OPD_R= sum(R_n, 3); 
refrence_phase= OPD_R*dX*2*pi/lambda;

%Sample beam
%Sample beam amplitude
S= ones(size_x, size_y); %amplitude of sample beam. assumed to be one and is irrelevant since angle
%Sample beam phase
%Preparing dimensions for axes

p_x=-size_x/2:(size_x/2-1);
p_y=-size_y/2:(size_y/2-1);

[x,  ~]=meshgrid(p_x*pix_ccd,p_y*pix_ccd);  % 2D space. wil be used later

%off axis: from B. A. E. Saleh and M. C. Teich, “Fourier optics,” in 
%Fundamentals of Photonics, B. A. E. Saleh ed. 2, pp. 140-141
%object wave S(x,y) is doubled by exp(-j*k*x*sin(alpha)) when the
%propagation is in the z plane. 
%The multiplication in this exponent in the x coordiante, unlike exp(j*fi)
%which is a constant, makes the cross- correlation terms be convoloved in
%the 2D x-y FT domain with the corresponding delta function, causing its
%deflection in the x axis.  
%U(x,y)=abs(R(x,y)+S(x,y)*exp(-j*k*x*sin(alpha))) %SIND !!!!! not SIN
off_axis_phase= (2*pi/lambda)*sind(alpha).*x;
%quadratic_phase= (x.^2+y.^2)/z; %a product of beam curvature- not a perfect plane wave

%For rotation
OPD_teta_array= zeros(size_x, size_y, N); %cell 1 is for 0 deg,.. cell N for 175 deg. 
for t=1:size_x %For each slice we rotate seperatly and than we sum the results to the
%corrsponding line in the OPD for each angle.
for l=1:N %For each angle
rot=  imrotate((squeeze(n(t,:,:))'), teta(l) , 'nearest','crop'); 
rot(rot==0)=medium;
OPD_teta_array(t,:, l)= sum(rot );
end
end

%Interfernce on the detector with no sample, but rather with a cuvatte
%filled with medium
n_cuv= medium*ones(size(n)); %goes through medium instead of sample. 
OPD_cuv= sum(n_cuv, 3); 
cuv_phase= OPD_cuv*dX*2*pi/lambda;

for k=1:N %For each angle
sample_phase_teta= OPD_teta_array(:,:,k)*dX*2*pi/lambda; %kz

%Interfernce on the detector with sample rotated in angle teta
%I(:,:,k+1)= (abs(R.*exp(1i*refrence_phase) + S.*exp(1i*quadratic_phase).*exp(1i*off_axis_phase).*exp(1i*sample_phase_teta))).^2;
I(:,:,k+1)= (abs(R.*exp(1i*refrence_phase) + S.*exp(1i*off_axis_phase).*exp(1i*sample_phase_teta))).^2;
end

%I(:,:,1)= (abs(R.*exp(1i*refrence_phase) + S.*exp(1i*quadratic_phase).*exp(1i*off_axis_phase).*exp(1i*cuv_phase))).^2;
I(:,:,1)= (abs(R.*exp(1i*refrence_phase) + S.*exp(1i*off_axis_phase).*exp(1i*cuv_phase))).^2;

end


