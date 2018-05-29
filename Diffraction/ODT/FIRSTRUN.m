%For first run with data
%Must pad with zeros to get back to normal size of picture (A algorithm)
%the output [x_c, y_c] result is an input for COMP
clear all
close all
clc

%%
%Optical set-up constants& parameters
% lambda= 632.8*10^(-9); %wavelength in [m]
lambda= 632.8*10^(-9);
SIZE= 256;
x_ap_c = 712; %Initial approximation
y_ap_c = 598;

%%
%OPDs extraction
%note that all multliplication in constants such as dx in the part prior
%to the beam refrencing is redundant and can be avoided. 
temp2 = imread('C:\Users\Gili\Documents\MATLAB\Yonina\T_cell\T_xy_1.png');
I(:,:,1)= temp2(y_ap_c-SIZE/2 :y_ap_c+ SIZE/2-1 ,x_ap_c-SIZE/2 :x_ap_c+ SIZE/2-1); %First cut

for i=1:73
temp2 = imread(strcat('C:\Users\Gili\Documents\MATLAB\Yonina\T_cell\', num2str(i),'.png'));
% temp2= rgb2gray(temp);
I(:,:,i+1)= temp2(y_ap_c-SIZE/2 :y_ap_c+ SIZE/2-1 ,x_ap_c-SIZE/2 :x_ap_c+ SIZE/2-1); %First cut
end

[rows,cols]= size(I(:,:,1));
N=cols/4; 
M=rows/4;
I_1= zeros(rows, cols,  size(I,3));
I_2= zeros(M, N,  size(I,3));
I_3= zeros(rows, cols,  size(I,3));
A= zeros(rows, cols,  size(I,3)-1);
x_c= zeros(size(I,3)-1,1);
y_c= zeros(size(I,3)-1,1);

for k=1: size(I,3) 
I_1(:,:,k)= fftshift(fft2(ifftshift(I(:,:,k)))); %fftshift vs ifftshift only matters for odd grids

%cc cropping
if (k==1)
    cc_location = find_max(I_1(SIZE/2+1,:,k) , cols);
end

I_2(:,:,k) = I_1((1.5*N+1):2.5*N, cc_location-rows/8:cc_location+rows/8 -1 ,k);

%%zero padding
I_3(1.5*N+1:2.5*N,1.5*N+1:2.5*N ,k)= I_2(:,:,k);

I_3(:,:,k) = fftshift(ifft2(ifftshift(I_3(:,:,k))));
 
I_3(:,:,k)= atan2 (  imag(I_3(:,:,k))  , real(I_3(:,:,k))  );

end

for k=1: size(I,3) -1
A(:,:,k)= I_3(:,:,k+1) - I_3(:,:,1);
A(:,:,k)=phase_unwrapping_2D(A(:,:,k));
end  


A=lambda*A/(2*pi);%Sign depends on positive/negative off axis angle
A(A<0)=0;
figure
imagesc (A(:,:,1))

for k=1: size(I,3) -1
[col_cm, row_cm]= find_COM(A(:,:,k)); %Tells us where to crop
x_c(k)= x_ap_c + col_cm-(SIZE/2 +1); %NEW center
y_c(k)= y_ap_c + row_cm-(SIZE/2 +1); %NEW center
end

for i=1:36
figure
imagesc(A(:,:,i))
end
