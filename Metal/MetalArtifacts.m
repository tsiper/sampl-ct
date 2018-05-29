%% Created metal artifact of Computed Tomography   %%
% This simulation is parallel beam CT.
% For metal artifacts reduction
% Titipong kaewlek M.Sc.(Radiological Sciences)
% Department of Radiological Technology                  
% Faculty of Allied Health Sciences 
% Phitsanulok 65000, Thailand ,e-mail:titipongk@nu.ac.th  


close all
clear
clc

e= [1   .920  .690   	0     	0        90
    -.8 	.874 	.6624   	0      -.0184     90
    -.2 	.310 	.1100     .22    	  0       72
    -.2 	.410 	.1600    -.22    	  0      108 
    .1  	.250 	.2100   	0     	.35       90 
    .1  	.0460 	.046   	    0     	.1         0 
    .1  	.0460 	.046  	    0      -.1         0 
    .1  	.0460 	.023     -.08      -.605       0 
    .1  	.0230 	.023    	0      -.605       0 
    .1  	.046  	.023   	  .06      -.605      90 

 
 %{ %} 
    % for small circle

    10  	.018 	.018   	0     	.35         0 % metal inside write region
    10      .018 	.018    -.22    0         90 % metal inside back region
    10      .018 	.018     .22    0         90 % metal inside back region
    10  	.018  	.018     0     -.35          0]; % metal on bottom 
 
    
%{
 % for big circle
    10  	.05 	.05   	0     	.35         0 % metal inside write region
    10      .05 	.05    -.22    0         90 % metal inside back region
    10      .05 	.05     .22    0         90 % metal inside back region
    10  	.05  	.05     0     -.35          0]; % metal on bottom 
%}

%{
 % for small oval
    10  	.056 	.023   	0     	.35         180 % metal inside write region
    10      .056 	.023    -.22    0         90 % metal inside back region
    10      .056 	.023     .22    0         90 % metal inside back region
    10  	.056  	.023     0     -.45          180]; % metal on bottom 
%}
%{
 % for big oval
    10  	.076 	.023   	0     	.35         180 % metal inside write region
    10      .076 	.023    -.22    0         90 % metal inside back region
    10      .076 	.023     .22    0         90 % metal inside back region
    10  	.076  	.023     0     -.45          180]; % metal on bottom 
  %}  
%{
 % for big oval
 
    10      .076 	.026    -.22    0         90 % metal inside back region
    10  	.076  	.026     0     -.75          180]; % metal on bottom 
 %}
    
theta1 =0:1:179;

Pori = phantom(e,256);

figure,imshow(Pori),title('original with metal non recon');%('original with metal non recon')

%theta1 =0:1:359;%(180 projection) if used different in 2 degree ; 0:2:178(90 projection)

[R,xp]=radon(Pori,theta1);

figure,imshow(R,[],'Xdata',theta1,'Ydata',xp,...
            'InitialMagnification','fit')
xlabel('\theta (degrees)')
ylabel('x''')
colormap(gray),title('shepp logan original')



XD=double(Pori);
[mm nn]=size(XD);
I55=[];
I66=[];
% thresholding for segmentation metal
for ii=1:mm
for jj=1:nn
    if XD(ii,jj)>1.3 
        I55(ii,jj)=1;
    else
        I55(ii,jj)=0;
    end
    if XD(ii,jj)<1.3  
        I66(ii,jj)=1;
    else
        I66(ii,jj)=0;
    end


end
end
I5=I55.*XD;
I6=I66.*XD;
[R2,xp2]=radon(I5,theta1);
[R3,xp3]=radon(I6,theta1);

figure,imshow(I5),title('metal original');

[mm nn]=size(R2);


Rn=R2;
for iii=1:mm
for jjj=1:nn
    
    if R2(iii,jjj)>0
        Rn(iii,jjj)=1;
    else    
        Rn(iii,jjj)=0;
    end


end
end


[uu vv]=size(R2);

Rnew=R2;

for ii=1:uu
for jj=1:vv
    if R2(ii,jj)>0
        Rnew(ii,jj)=0;
    else
        Rnew(ii,jj)=1;
    end


end
end



Rm=max(max(R));



sinoPre=R(1:uu,1:vv);    

sinoNew = sinoPre;
R2n=[];
%{%}
%normalized data  on projection
for kk = 1:vv
    y1 = find(Rnew(:,kk)==0);
    n_y = length(y1);
    if n_y >1
        yy1 = 1:uu;

        R2n(y1,kk)=(50/100)*Rm;
        sinoNew(y1,kk) =R2n(y1,kk);
    
      
    end;
end;
%Image genarate artifacts

PNew= iradon(sinoNew, theta1,'linear','shepp-logan',0.9);

figure,imshow(PNew),title(' metal artifacts reconstruction');


%Copyright (c) 2012, Titipong Kaewlek
%All rights reserved.


%% NOGA & DOR

Artifacts=PNew(2:257,2:257)-Pori;
Artifacts_synogram=radon(Artifacts,theta1);
figure;
subplot(231);
imshow(Pori);
colormap('bone');
title('Original');
subplot(232);
imshow(PNew);
colormap('bone');
title('Artifacts');
subplot(233)
imshow(Artifacts);
colormap('bone');
title('Only artifacts');
subplot(234);
imagesc(sinoPre);
colormap('bone');
title('Original sinogram');
subplot(235);
imagesc(sinoNew);
colormap('bone');
title('Artifacts sinogram');
subplot(236);
imagesc(Artifacts_synogram);
colormap('bone');
title('Only artifacts synogram');

figure;
subplot(211);
imagesc(Artifacts_synogram);
colormap('bone');
title('Only artifacts synogram');
subplot(212);
imagesc(sinoNew-sinoPre);
colormap('bone');
title('Only artifacts synogram (reduction)');

figure;
subplot(211);
imshow(Pori);
colormap('bone');
title('Original');
subplot(212);
imshow(iradon(sinoNew-Artifacts_synogram,theta1));
colormap('bone');
title('Artifact Reduction');