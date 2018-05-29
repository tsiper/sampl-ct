function [ DecCoordinates,DecSinogram ] = DecimateSPURS( Coordinates, Sinogram, DecFactor )
%DECIMATESPURS This function decimates the sinogram by uniform angles according
%to the provided DecFactor

[m,n] = size(Sinogram);
DecSinogram = zeros(m,floor(n/DecFactor));
DecCoordinates = zeros( m*floor(n/DecFactor) , 2);

% Decimating the sinogram and coordinates
i = 1;
for j=1:DecFactor:n
    DecSinogram(:,i) = Sinogram(:,j);
    DecCoordinates((i-1)*m+1:i*m,:) = Coordinates((j-1)*m+1:j*m,:);
    i = i+1;
end

%


