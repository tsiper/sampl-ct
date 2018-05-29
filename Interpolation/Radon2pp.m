function [ ppSinogram ] = Radon2pp( Sinogram , thetas )
%RADON2PPRADON This function takes the Radon transform Sinogram and transforms
%it into a Pseudo-Polar Sinogram we assume thetas are equi-spaced

[m,~] = size(Sinogram);
n = (m-1)/2;


% We first perform an upsampling step, using prior information about the signal

% Our original polar grid points
x = -n:n;
thetas_rad = thetas / 180 * pi;
ppSinogram = zeros(size(Sinogram));

% Sinogram = F1(Sinogram);


for i=1:length(thetas_rad)
        % Our new grid, according to theta angle
        if (thetas_rad(i) >= -pi/4) && (thetas_rad(i) <= pi/4)
            xx = (x)*cos((thetas_rad(i)));
            ppSinogram(:,i) = spline(x,(Sinogram(:,i)),xx)*abs(cos(thetas_rad(i)));
%             ppSinogram(:,i) = interp1(x,Sinogram(:,i),xx,'pchip')*abs(cos(thetas_rad(i)));
        else
            xx = (x)*sin((thetas_rad(i)))-4*(thetas_rad(i)-pi/4)/pi;
            ppSinogram(:,i) = spline(x,(Sinogram(:,i)),xx)*abs(sin(thetas_rad(i)));
%             ppSinogram(:,i) = interp1(x,Sinogram(:,i),xx,'pchip')*abs(sin(thetas_rad(i)));
        end
end

% ppSinogram = real(invF1(ppSinogram));

end