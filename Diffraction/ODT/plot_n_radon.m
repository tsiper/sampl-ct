function [] = plot_n_radon(nr1,nr2,range,SIZE)

% Plotts 2 n_radon ouputs and compares them

figure
subplot(2,3,1)
imagesc(squeeze(nr1(:,SIZE/2+1,:)),range)
title 'sagittal slice- radon - Gili'
colorbar
hold on
subplot(2,3,2)
imagesc(squeeze(nr1(SIZE/2+1,:,:)),range)
title 'axial slice- radon - Gili'
colorbar
hold on
subplot(2,3,3)
imagesc(flipud(nr1(:,:,SIZE/2+1)),range)
title 'coronal slice- radon - Gili'
colorbar
hold on
subplot(2,3,4)
imagesc(squeeze(nr2(:,SIZE/2+1,:)),range)
title 'sagittal slice- radon - AT'
colorbar
hold on
subplot(2,3,5)
imagesc(squeeze(nr2(SIZE/2+1,:,:)),range)
title 'axial slice- radon - AT'
colorbar
hold on
subplot(2,3,6)
imagesc(flipud(nr2(:,:,SIZE/2+1)),range)
title 'coronal slice- radon - AT'
colorbar
hold off

%CompareImages(nr1,nr2);

end

