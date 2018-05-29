function [n2, n1] = tomo_rytov( u_rytov, lambda, dX, theta, nm, iter_num)

%%PARAMETERS
[~,N,L]= size(u_rytov);
%Beam parameters
k0=2*pi/lambda;
%Optical system parameters
dk=(2*pi/N)*(1/dX);

%%
%%INITIALIZATION
F=single(zeros(N,N,N)); %Fourier space of f
count=single(zeros(N,N,N)); 

%%
%%MAPPING THE SCATTERED WAVES
%Incident plane wave
ki_x=0; ki_y=0;
ki_z=nm*k0;

%Scattered waves
P=(1:N)-(N/2+1);
[ks_x,ks_y]=meshgrid(P*dk,P*dk);
ks_z=real(sqrt((nm*k0)^2-ks_x.^2-ks_y.^2)); %ki_x.^2+ ki_y^2+ki_z^2=ks_x.^2+ ks_y^2+ks_z^2

%Eliminate unnecessary waves (evanescent waves)
[xx,yy]=find((ks_z)>0);
ks_x=ks_x(min(xx):max(xx),min(yy):max(yy)); 
ks_y=ks_y(min(xx):max(xx),min(yy):max(yy)); 
ks_z=ks_z(min(xx):max(xx),min(yy):max(yy)); 
ks_zc=ks_z; %copy of ks_z
ks_x=ks_x(:)'; ks_y=ks_y(:)'; ks_z=ks_z(:)';

%% Extracting Phase maps

%%Finding the rytov field and MAPPING ON EWALD SPHERE
for jj=1:L
    U_rytov=dX^2.*fftshift(fft2(ifftshift(u_rytov(:,:,jj))));  
    F_2D=(1i*2*ks_zc).*U_rytov(min(xx):max(xx),min(yy):max(yy)); %(1i*ks_zc/pi).*U_rytov(min(xx):max(xx),min(yy):max(yy));

    %Do the rotation
    R=[cosd(-theta(jj)) -sind(-theta(jj)); sind(-theta(jj)) cosd(-theta(jj))]; %Rotation matrix that rotates the scattered and incident waves according to theta    
    ks_r=R*[(ks_x-ki_x); (ks_z-ki_z)];
    ks_xr=ks_r(1,:); ks_zr=ks_r(2,:);
    
    % rotating around y-> x and z change.
    % Choice of order (rows/cols) doesn't matter since the only differnce is
    % +/-sin(teta) and teta goes from -pi to pi so it is symmetric
    % --> equivelent to:    
    %     ks_r=R*[(ks_z-ki_z); (ks_x-ki_x)]; 
    %     ks_xr=ks_r(2,:); ks_zr=ks_r(1,:);

    %mapping the rotated scattered waves to indices
    indKy=round((ks_y)/dk+(N/2+1)); 
    indKz=round((ks_zr)/dk+(N/2+1));
    indKx=round((ks_xr)/dk+(N/2+1));
    F_2D=F_2D(:);
    
    h=find(indKz>N | indKx>N | indKz<1 | indKx<1);
    indKx(h)=[]; indKy(h)=[]; indKz(h)=[]; F_2D(h)=[];
    
    ind3D=indKy+(indKx-1)*N+(indKz-1)*N^2;
    F(ind3D)=F(ind3D)+F_2D';
    count(ind3D)=count(ind3D)+1;
          
end
clear u_rytov phase complex_field ks_x ks_y ks_z ks_zc ks_yr ks_zr indKx indKy indKz F_2D ind3D U_rytov ks_r ks_xr xx yy
%%
%%FINAL PROCESSING
i=find(count>0);
F(i)=F(i)./count(i);
clear count i 

f=fftshift(ifftn(ifftshift(F)))*(N*dk/(2*pi))^3;
n=(sqrt(nm^2-f./((k0).^2)));
n1= real(n);
n2=n1;

%Iterative non-negitivity constraint alghorithm- corrected
for i=1:iter_num
n2(n2<nm)=nm;
f=(k0^2)*(nm^2-n2.^2);
tempF=fftshift(fftn(ifftshift(f)))*(1/(N*dk/(2*pi))^3); % new 3D Fourier field
tempF(abs(F)~=0)= F(abs(F)~=0);
f=fftshift(ifftn(ifftshift(tempF)))*(N*dk/(2*pi))^3;
n2=real(sqrt(nm^2-f./((k0).^2)));
end

end
