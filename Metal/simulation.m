
clear all;
% close all;
clc;

%% Create Phantom:

res = 256;
phantom = zeros(res);
numOfMaterials = 4;

x = repmat(1:res,res,1);
y = repmat((1:res)',1,res);

circleCenterX = [120,90,170,100];
circleCenterY = [120, 85, 125,150];
circleRadius = [90, 20, 15,25];

for i=1:numOfMaterials
    z = (x-circleCenterX(i)).^2 + (y-circleCenterY(i)).^2;
    phantom(z <= circleRadius(i)^2) = i;
end

figure;
imagesc(phantom);

%% Choose Materials:

materials = {'soft_tissue','Ti','Fe','blood'};


%% Load spectra:

load 'spectra.mat';

%% Attenuations Coeffs:

attCoeffs = cell(numOfMaterials,1);

figure;
for i=1:numOfMaterials
   attCoeffs{i} = getMassAttenCoeff(materials{i},Spectra{1}.energies);
   subplot(1,numOfMaterials,i);
   plot(Spectra{1}.energies,attCoeffs{i});
   title(materials{i});
end

%% Scan:
d_theta = 1;
theta = 1:d_theta:180;
PixSize = 0.00001; % [m]
EngSize = sum(diff(Spectra{1}.energies));
sinogram = zeros(size(theta,1),res);

figure;
for i=1:length(theta)
   pic = imrotate(phantom,theta(i));
   l = size(pic,1);
   d = ceil((l - res)/2);
   if mod(l,2)
        phantomForScan = pic(d+1:end-d+1,d+1:end-d+1);
   else
        phantomForScan = pic(d+1:end-d,d+1:end-d);
   end
   
%    imagesc(phantomForScan);
%    title([num2str(i) ',' num2str(size(phantomForScan,1))]);
%    pause(0.00001);
   
   for m = 1:res
      current_spectra = Spectra{1}.Intensity';
      for n = 1:res
          if (phantomForScan(m,n))
                current_spectra = current_spectra.*exp(-sum(attCoeffs{phantomForScan(m,n)})*PixSize);
          end
      end
      sinogram(i,m) = EngSize * sum(current_spectra);
   end
end
sinogram =sinogram';


figure;
imagesc(sinogram);

%% Adjustmens:

adjusted_sinogram = -(sinogram - max(max(sinogram)));
reconstructedPhantom = iradon((adjusted_sinogram),theta);
% reconstructedPhantom = iradon(-log(sinogram),theta);

figure;
imagesc(reconstructedPhantom);
