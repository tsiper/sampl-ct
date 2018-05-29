clear
load('run_6.mat');
runInd = 2;
titles = {'run 1: more circles with noise',...
        'run 2: circles with noise',...
        'run 3: 2 bins 3 materials without noise',...
        'run 4: water and adipose with noise (overlapping)',...
        'run 5: I and Ca with noise (overlapping)',...
        'run 6: 3 materials where only 2 are reconstrcuted',...
        'run 7: high resolution with noise level of 5%',...
        'run 8: reconstruction from a switching scan',...
        };
%% plot phantom
if I == 3
    myRGBphantom = RGBphantom;
elseif I == 2
    myRGBphantom = cat(3,RGBphantom, zeros(length(RGBphantom(:,1,1))));
end
% figure; imshow(myRGBphantom);
% title(chosen_phantom);
%% plot projections with noise
if k == 3
    RGBprojections = cat(3,p_input{1},p_input{2},p_input{3});
elseif k == 2
    RGBprojections = cat(3,p_input{1},p_input{2},zeros(size(p_input{1})));
end
RGBprojections = (RGBprojections-min(RGBprojections(:)))/range(RGBprojections(:));
figure; imshow(RGBprojections);
title(['Projections, scantype is ',scan_type]);

%% plot reconstructed image
if I == 3
    RecRGBmaterials = cat(3,Hpp_T(a_hat{1})/ReconstructionParams.base(1),...
        Hpp_T(a_hat{2})/ReconstructionParams.base(2),...
        Hpp_T(a_hat{3})/ReconstructionParams.base(3));
elseif I == 2
    RecRGBmaterials = cat(3,Hpp_T(a_hat{1})/ReconstructionParams.base(1),...
    Hpp_T(a_hat{2})/ReconstructionParams.base(2),...
    zeros(size(Hpp_T(a_hat{1}))));
elseif IRecon == 2
    RecRGBmaterials = cat(3,zeros(size(Hpp_T(a_hat{1}))),...
        Hpp_T(a_hat{1})/ReconstructionParams.base(1),...
    Hpp_T(a_hat{2})/ReconstructionParams.base(2));
end
figure;
subplot(131); imshow(myRGBphantom); title('original phantom');
subplot(132); imshow(RecRGBmaterials); title('reconstructed phantom');
subplot(133); imshow(myRGBphantom-RecRGBmaterials,[]); title('difference');
suptitle(titles{runInd});
%% ploting simulator results:
figure;

subtightplot(321); imshow(a_hat{1},[]); title('Material Sinogram: Iodine');
subtightplot(322); imshow(a_hat{2},[]); title('Material Sinogram: Calcium');

subtightplot(323); imshow(Hpp_T(a_hat{1}),[]); title('Material Map: Iodine');
subtightplot(324); imshow(Hpp_T(a_hat{2}),[]); title('Material Map: Calcium');

subtightplot(325); imagesc(phantom{1}-Hpp_T(a_hat{1})); title('First Material Difference: Iodine'); axis off; colorbar('EastOutside');
subtightplot(326); imagesc(phantom{2}-Hpp_T(a_hat{2})); title('Second Material Difference: Calcium');axis off; colorbar('EastOutside');

%% Fista Run
openfig('RLcircles_lowRes_FISTA.fig');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear
% low resolution results of RLcircles
load('zubal64_data.mat');

%% plot phantom, projactions, and initial reconstructions
figure('Position',[50 50 1400 800]);

RGBradonMaterials = cat(3,3*RadonMaterials{1},0.4*RadonMaterials{2},0.2*RadonMaterials{3});
range_radonMaterials = range(RGBradonMaterials(:));
min_radonMaterials = min(RGBradonMaterials(:));
RGBradonMaterials = (RGBradonMaterials-min_radonMaterials)/range_radonMaterials;
RGBmuEnergies = cat(3,Hpp_T(p_input{1}),Hpp_T(p_input{2}),Hpp_T(p_input{3}));
RGBmuEnergies = (RGBmuEnergies-min(RGBmuEnergies(:)))/range(RGBmuEnergies(:));
base = ReconstructionParams.base;
RecRGBradonMaterials = cat(3,3*a_hat{1},0.4*a_hat{2},0.2*a_hat{3});
RecRGBradonMaterials = (RecRGBradonMaterials-min_radonMaterials)/range_radonMaterials;
RecRGBmaterials =cat(3,...
    Hpp_T(a_hat{1}/base(1)),...
    Hpp_T(a_hat{2}/base(2)),...
    Hpp_T(a_hat{3}/base(3))...
    );

% Preperaing for figure- Shahar
RGBradonMaterials    = trimcols(RGBradonMaterials,floor(size(RGBradonMaterials,1)/2));
RecRGBradonMaterials = trimcols(RecRGBradonMaterials,floor(size(RecRGBradonMaterials,1)/2));
RGBmuEnergies        = rgb2gray(RGBmuEnergies);
GraySinogram         = trimcols(rgb2gray(cat(3,p_input{1},p_input{2},p_input{3})),floor(size(p_input{1},1)/2));


% ax(1) = subtightplot(231); imshow(RGBradonMaterials); title('Materials Sinograms');
ax(2) = subtightplot(232); imshow(GraySinogram); title(sprintf('Thresholding:\n50, 80, 140 keV'));

ax(3) = subtightplot(233); imshow(RecRGBradonMaterials); title('Reconstructed Sinograms');

ax(4) = subtightplot(2,3,[1 4]); imshow(RGBphantom); title(sprintf('Materials Phantom:\nI, Ca, Soft Tissue'));
ax(5) = subtightplot(235); imshow(RGBmuEnergies); title('Energy Attenuation Map');

ax(6) = subtightplot(236); imshow(RecRGBmaterials); title('Reconstructed Phantom');

set(ax,'FontSize',22);

openfig('zubal64_FISTA.fig');

SaveFigure('Spectral_Zubal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear
% low resolution results of RLcircles
load('RLcirc64_data.mat');

%% plot phantom, projactions, and initial reconstructions
figure('Position',[50 50 1400 800]);

RGBradonMaterials = cat(3,3*RadonMaterials{1},0.4*RadonMaterials{2},0.2*zeros(size(RadonMaterials{1})));
range_radonMaterials = range(RGBradonMaterials(:));
min_radonMaterials = min(RGBradonMaterials(:));
RGBradonMaterials = (RGBradonMaterials-min_radonMaterials)/range_radonMaterials;
RGBmuEnergies = cat(3,Hpp_T(p_input{1}),Hpp_T(p_input{2}),Hpp_T(zeros(size(p_input{2}))));
RGBmuEnergies = (RGBmuEnergies-min(RGBmuEnergies(:)))/range(RGBmuEnergies(:));
base = ReconstructionParams.base;
RecRGBradonMaterials = cat(3,3*a_hat{1},0.4*a_hat{2},zeros(size(a_hat{1})));
RecRGBradonMaterials = (RecRGBradonMaterials-min_radonMaterials)/range_radonMaterials;
RecRGBmaterials = cat(3,Hpp_T(a_hat{1}/base(1)),Hpp_T(a_hat{2}/base(2)),Hpp_T(zeros(size(a_hat{2}))));

subplot(231); imshow(cat(3,p_input{1},p_input{2},zeros(size(p_input{2})))); title(sprintf('Energy Sinograms:\n70, 140 keV'));
subplot(232); imshow(RGBradonMaterials); title('Materials Sinograms');
subplot(233); imshow(RecRGBradonMaterials); title('Reconstructed Sinograms');

subplot(234); imshow(RGBmuEnergies); title('Effective Attenuation of the Energies');
subplot(235); imshow(cat(3,RGBphantom,zeros(size(RGBphantom(:,:,1))))); title(sprintf('Materials Phantom:\nI, Ca'));
subplot(236); imshow(RecRGBmaterials); title('Reconstructed Phantom');

figure; imshow(cat(3,RGBphantom,zeros(size(RGBphantom(:,:,1))))-RecRGBmaterials,[]);
% openfig('RLcircles_lowRes_FISTA.fig');

