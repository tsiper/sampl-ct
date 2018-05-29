
N = 500;
proj_num = 400;
thetas = linspace(0,2*pi,proj_num);
d_theta = thetas(2)-thetas(1);
P = 500;
%% load projections
% raw_data = zeros(N,N,proj_num);
f = figure;
tic
for n = 0:proj_num-1
%     if n<10
%         n_str = ['00',num2str(n)];
%     elseif n<100
%         n_str = ['0',num2str(n)];
%     else
%         n_str = num2str(n);
%     end
%     path_str = ['projections_35\projection0',n_str,'.tif'];
%     raw_data(:,:,n+1) = imread(path_str);
    im = raw_data_35(:,:,n+1);
    figure(f);
    imshow(im,[]);
    drawnow;
    pause(0.01);
end
%%
toc
min_val = min(raw_data_35(:));
norm_data =1- (raw_data_35-min_val)/range(raw_data_35(:));
rec = zeros(280,280,N);
r = figure();
display('Performs Reconstruction');
tic
% fanbeam
for m = 1:N
    
    F = shiftdim(raw_data_35(m,1:end,1:end),1);
%     figure; 
%     subplot(221);imagesc(F);
%     R = fan2para(F,30000,'FanCoverage','cycle','FanRotationIncrement',d_theta*180/pi,...
%          'FanSensorGeometry','line','FanSensorSpacing',1, 'ParallelCoverage','cycle',...
%          'ParallelRotationIncrement',d_theta*180/pi,'ParallelSensorSpacing',1);
%      subplot(222);imagesc(R);
%      t = size(R,2)
%     rec2 = iradon(R(:,1:round(t/2)),round(t/2));
%     rec1 = ifanbeam(F,30000,'FanCoverage','cycle','FanRotationIncrement',d_theta*180/pi,...
%          'FanSensorGeometry','line','FanSensorSpacing',1);
% %     rec(:,:,m) = rec1;
% %     imwrite(255*(rec(:,:,m)),['C:\Users\student\Downloads\phantom\rec\rec',num2str(m),'.tiff']);
%     subplot(223);
%     imshow(rec1,[]);
%     subplot(224);
%     imshow(rec2,[]);
%     drawnow;
d=30000;
   [P,oploc,optheta] = fan2para(F,d,...
                             'FanSensorSpacing',1,...
                             'ParallelSensorSpacing',1,...
                             'FanSensorGeometry','line',...
                             'Interpolation','spline',...
                             'FanCoverage','cycle',...
                             'FanRotationIncrement',d_theta*180/pi); 

optional_args = {'spline', 'Ram-Lak'};
[I,H] = iradon(P,optheta);%,optional_args{:});
figure(r);imagesc(I)
% display(m)
drawnow;
end
% p_35.p = P
% p_35.thetas = optheta
% save('radonProjection.mat','p_20','p_28','p_35');
%%
% threshold = 0.00438;
% threshold = 0.0043;
% 
% bin_rec = zeros(size(rec));
% bin_rec(rec>=threshold) = 1;
% points = find(bin_rec==1);
% %%
% vol = linspace(0,1,N);
% [X,Y,Z] = meshgrid(vol,vol,vol);
% figure;
% scatter3(X(points),Y(points),Z(points),'.')
camera_size = 0.05;
pix_num = 1024;
pix_size = camera_size/pix_num;
rot_distance = 39.5-28;
dot_dist_in_pix = rot_distance/pix_size
toc
%%
n=30;
filt = 1/n*ones(n,1);
out = conv(GilW350dot01mA.Impulses,filt);
plot(GilW350dot01mA.Channelnumber,GilW350dot01mA.Impulses)
hold on
plot(GilW350dot01mA.Channelnumber,out(n/2:length(GilW350dot01mA.Channelnumber)+n/2-1));
%%
norm_spec_20 = spec_20/sum(spec_20);
norm_spec_28 = spec_28/sum(spec_28);
norm_spec_35 = spec_35/sum(spec_35);
figure;
plot(1:length(norm_spec_20),mean(mean(mean(raw_data_20(1:40,1:40,:))))*norm_spec_20,...
    1:length(norm_spec_28),mean(mean(mean(raw_data_28(1:40,1:40,:))))*norm_spec_28,...
    1:length(norm_spec_35),mean(mean(mean(raw_data_35(1:40,1:40,:))))*norm_spec_35);

save('spectras.mat','spec_20','spec_28','spec_35');