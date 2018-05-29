function [ RBC ] = createRBC_disocyte( size, dX)
%based on Formula #3 in Zeidan Adel, and Dvir Yelin. "Reflectance confocal microscopy of
%red blood cells: simulation and experiment." Biomedical optics express 6.11 (2015): 4335-4343., formula 3

%Tonicity [mO]- a measure of the effective osmotic pressure gradient

%coefficients 
R0 = 4*10^-6; %[m] 
C0 = 0.35*10^-6; %[m]
C2 = 3.8*10^-6; %[m]
C4 = -3.5*10^-6; %[m]

p_x=-size/2:(size/2-1);
p_y=-size/2:(size/2-1);
[X,  Y]=meshgrid(p_x*dX,p_y*dX);

height = real((  C0  +  C2*((X.^2+Y.^2)/R0^2) + C4*((X.^2+Y.^2).^2/R0^4)   ) .* sqrt(1 - ((X.^2+Y.^2)/R0^2)));

RBC= zeros(size,size,size);

for i=1:size
    for j=1:size
        L= round(height(i,j)/dX);
        RBC(i,j,size/2-(L-1):size/2+1+L-1)=1;     
    end
end


end

