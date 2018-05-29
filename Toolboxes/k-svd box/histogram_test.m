close all

temp_scanlow=dicomread('sheepscanlow/sheepscanlow00001.dcm');
temp_scanlow=double(temp_scanlow);
maxim=max(max(temp_scanlow));
temp_scanlow=temp_scanlow./maxim;
hist=histogram(temp_scanlow,100);
adjusted_image=imadjust(temp_scanlow,[0.15,0.37],[0,0.25]);
%figure();imshow(temp_scanlow,[]);
figure();imshow(adjusted_image,[]);

temp_scanhigh=dicomread('sheepscanhigh/sheepscanhigh00001.dcm');
temp_scanhigh=double(temp_scanhigh);
maxim=max(max(temp_scanhigh));
temp_scanhigh=temp_scanhigh./maxim;
figure();
hist_high=histogram(temp_scanhigh,100);
adjusted_image_high=imadjust(temp_scanhigh,[0.15,0.37],[0,0.25]);
%figure();imshow(temp_scanhigh,[]);
figure();imshow(adjusted_image_high,[]);