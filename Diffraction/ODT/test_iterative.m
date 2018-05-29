clear all
close all
clc

%% experimetal set-up
pix_ccd= 5.2e-06; %[m];
Mag= 75;
dX= pix_ccd/Mag;
nm= 1.3334; %medium
lambda=633*10^(-9); %wavelength in [m]- probably, maybe 633
SIZE= 256;
K=2*pi/lambda;

% simulation data creation
m = SIZE;
n = SIZE;
teta= -60:60;
k= length(teta);
nRBC=1.395;
n_diff= nRBC-nm;
RBC = n_diff*createRBC_disocyte(SIZE, dX); %for normal cells(discocytes)
proj_arr= rbcROTATE(RBC, teta, SIZE, dX, K );
IT_NUM=100;

%tomogrpahy- radon
    temp_RI = tomo_radon_2D(proj_arr,lambda, dX,teta, nm );     
    RI_NO_IT= temp_RI;
    teta_vec=-180:180;
    for it=1:IT_NUM %iterations
        %Revise unrealistic data
        temp_RI(temp_RI<nm)=nm;
        %Radon transform.
        new_arr = create_projs(temp_RI, lambda, dX, teta_vec, nm);
        %Replacement.
        for te= 1:length(teta)
        rep=  find(teta_vec==teta(te));
        new_arr(:,:,rep)= proj_arr(:,:,te);
        end 
        %Inverse Radon transform+Revise unrealistic data
        temp_RI = tomo_radon_2D(new_arr,lambda, dX,teta_vec, nm );  
    end
    RI_IT= temp_RI;
    
%tomogrpahy- rytov
[ODT_RI_IT,ODT_RI_NO_IT] = tomo_rytov(1i*proj_arr,lambda, dX,teta, nm, IT_NUM);     
    
RI=RBC+nm;
     
    figure
    subplot(3,1,1)
    imagesc(RI_NO_IT(:,:,SIZE/2+1),[1.33 1.395])
    colorbar
    hold on
    subplot(3,1,2)
    imagesc(RI_IT(:,:,SIZE/2+1),[1.33 1.395])
    colorbar
    hold on
    subplot(3,1,3)
    imagesc(RI(:,:,SIZE/2+1),[1.33 1.395])
    colorbar
    
    figure
    subplot(3,1,1)
    imagesc(squeeze(RI_NO_IT(:,SIZE/2+1,:)),[1.33 1.395])
    colorbar
    hold on
    subplot(3,1,2)
    imagesc(squeeze(RI_IT(:,SIZE/2+1,:)),[1.33 1.395])
    colorbar
    hold on
    subplot(3,1,3)
    imagesc(squeeze(RI(:,SIZE/2+1,:)),[1.33 1.395])
    colorbar
    
    figure
    subplot(3,1,1)
    imagesc(squeeze(RI_NO_IT(SIZE/2+1,:,:)),[1.33 1.395])
    colorbar
    hold on
    subplot(3,1,2)
    imagesc(squeeze(RI_IT(SIZE/2+1,:,:)),[1.33 1.395])
    colorbar
    hold on
    subplot(3,1,3)
    imagesc(squeeze(RI(SIZE/2+1,:,:)),[1.33 1.395])
    colorbar
        
    figure
    subplot(3,1,1)
    imagesc(ODT_RI_NO_IT(:,:,SIZE/2+1),[1.33 1.395])
    colorbar
    hold on
    subplot(3,1,2)
    imagesc(ODT_RI_IT(:,:,SIZE/2+1),[1.33 1.395])
    colorbar
    hold on
    subplot(3,1,3)
    imagesc(RI(:,:,SIZE/2+1),[1.33 1.395])
    colorbar
    
    figure
    subplot(3,1,1)
    imagesc(squeeze(ODT_RI_NO_IT(:,SIZE/2+1,:)),[1.33 1.395])
    colorbar
    hold on
    subplot(3,1,2)
    imagesc(squeeze(ODT_RI_IT(:,SIZE/2+1,:)),[1.33 1.395])
    colorbar
    hold on
    subplot(3,1,3)
    imagesc(squeeze(RI(:,SIZE/2+1,:)),[1.33 1.395])
    colorbar
    
    figure
    subplot(3,1,1)
    imagesc(squeeze(ODT_RI_NO_IT(SIZE/2+1,:,:)),[1.33 1.395])
    colorbar
    hold on
    subplot(3,1,2)
    imagesc(squeeze(ODT_RI_IT(SIZE/2+1,:,:)),[1.33 1.395])
    colorbar
    hold on
    subplot(3,1,3)
    imagesc(squeeze(RI(SIZE/2+1,:,:)),[1.33 1.395])
    colorbar
