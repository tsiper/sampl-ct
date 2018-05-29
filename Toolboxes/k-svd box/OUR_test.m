%  close all
%  clear all


%% read images %%
%%readind the original image%%
im_1 = dicomread('sheepscanlow/sheepscanlow00001.dcm');
im_1=double(im_1);
maximL=max(max(im_1));
im_1=im_1./maximL;

im_1=imcrop(im_1,[100 150 250 250]);
im_1=imadjust(im_1,[0.15,0.37],[0,0.25]);


%%creating the first image of the training set of high def%%
im_2 = dicomread('sheepscanhigh/sheepscanhigh00002.dcm');
im_2=double(im_2);
maximH=max(max(im_2));
im_2=im_2./maximH;

% im_2=imcrop(im_2,[100 150 250 250]);
im_2=imadjust(im_2,[0.15,0.37],[0,0.25]);

%image registration stage%
Tl=dicomread('sheepscanlow/sheepscanlow00003.dcm'); 
Tl=double(Tl);
Tl=Tl./max(max(Tl));


Tl=imadjust(Tl,[0.15,0.37],[0,0.25]);

Th=RegisterImages(Tl,im_2);
Tl=imcrop(Tl,[100 150 250 250]);
Th=imcrop(Th,[100 150 250 250]);




    


%%creating the training set%%
for i=3:7
        im_addrH=strcat('sheepscanhigh/sheepscanhigh0000', num2str(i),'.dcm');
        im_tempH = dicomread(im_addrH);
        im_tempH = double(im_tempH);
        maximtemp=max(max(im_tempH));
        im_tempH=im_tempH./maximtemp;

        im_addrL=strcat('sheepscanlow/sheepscanlow0000', num2str(i+2),'.dcm');
        im_tempL = dicomread(im_addrL);
        im_tempL = double(im_tempL);
        maximtemp=max(max(im_tempL));
        im_tempL=im_tempL./maximtemp;
        

        im_tempL=imadjust(im_tempL,[0.15,0.37],[0,0.25]);
        im_tempH=imadjust(im_tempH,[0.15,0.37],[0,0.25]);
        
        im_tempH=RegisterImages(im_tempL,im_tempH);
        
        im_tempL=imcrop(im_tempL,[100 150 250 250]);
        im_tempH=imcrop(im_tempH,[100 150 250 250]);
        
        Tl= [Tl, im_tempL];
        Th=[Th,im_tempH];
end

for i=10:25
 
        im_addrH=strcat('sheepscanhigh/sheepscanhigh000', num2str(i),'.dcm');
        im_tempH = dicomread(im_addrH);
        im_tempH = double(im_tempH);
        maximtemp=max(max(im_tempH));
        im_tempH=im_tempH./maximtemp;

        im_addrL=strcat('sheepscanlow/sheepscanlow000', num2str(i+2),'.dcm');
        im_tempL = dicomread(im_addrL);
        im_tempL = double(im_tempL);
        maximtemp=max(max(im_tempL));
        im_tempL=im_tempL./maximtemp;
        
        im_tempL=imadjust(im_tempL,[0.15,0.37],[0,0.25]);
        im_tempH=imadjust(im_tempH,[0.15,0.37],[0,0.25]);
        
        im_tempH=RegisterImages(im_tempL,im_tempH);
        
        im_tempL=imcrop(im_tempL,[100 150 250 250]);
        im_tempH=imcrop(im_tempH,[100 150 250 250]);
        
        Tl= [Tl, im_tempL];
        Th=[Th,im_tempH];
end

% for i=100:121
%  
%         im_addr=strcat('sheepscanhigh/sheepscanhigh00', num2str(i),'.dcm');
%         im_tempH = dicomread(im_addr);
%         im_tempH = double(im_tempH);
%         maximtemp=max(max(im_tempH));
%         im_tempH=im_tempH./maximtemp;        
% %         x=TVD_mm(im_tempH,300,30);
% %         im_tempH=reshape(x,512,512);
%         Th= [Th, im_tempH];
% end







%% set parameters %%
params.x = im_1;
params.data = Tl;
params.Tdata = 2;
params.blocksize = 8;
params.dictsize = 1024;
params.maxval = max(max(im_1));
params.trainnum = 64000;
params.iternum = 10;
params.memusage = 'high';
params.muthresh= 0.99;
%params.exact= 1;




%% create training data %%
data = params.data;
blocksize = params.blocksize;
trainnum = params.trainnum;
dictsize = params.dictsize;

p = ndims(data);

ids = cell(p,1);
if (p==1)
  ids{1} = reggrid(length(data)-blocksize+1, trainnum, 'eqdist');
else
  [ids{:}] = reggrid(size(data)-blocksize+1, trainnum, 'eqdist');
end
params.data = sampgrid(Tl,blocksize,ids{:});

%% high dataset
high_dataset = sampgrid(Th,blocksize,ids{:});
%%
params.data=[params.data ; high_dataset];

%% run k-svd training %%
D= ksvd(params,'');

%%%     

Imtcol=im2colstep(im_1,[params.blocksize, params.blocksize]);

%%seperation of dictionaries to Low & High%%
D_low=D(1:params.blocksize*params.blocksize,:);
D_low=normcols(D_low);
D_high=D(params.blocksize*params.blocksize+1:params.blocksize*params.blocksize*2,:);
D_high=normcols(D_high);

%%Sparse coding stage%%
SparseMat=omp(D_low'*Imtcol,D_low'*D_low,params.Tdata);

%%Reconstruction%%
Y_low=D_low*SparseMat;
Y_low= col2imstep(Y_low,[size(im_1,1) size(im_1,1)],[params.blocksize params.blocksize]);

Y_high=D_high*SparseMat; 
Y_high= col2imstep(Y_high,[size(im_1,1) size(im_1,1)],[params.blocksize params.blocksize]);

 params.dict=D_low;
% params.x = Y;

%%matching images parametes for the PSNR calculation%%

Y_low=mat2gray(Y_low);
Y_high=mat2gray(Y_high);

im_clean= dicomread('sheepscanhigh/sheepscanhigh00001.dcm');
im_clean=double(im_clean);
im_clean=im_clean./max(max(im_clean));

im_clean=imcrop(im_clean,[100 150 250 250]);
im_clean=imadjust(im_clean,[0.15,0.37],[0,0.25]);

psnr_Yh=psnr(Y_high,im_clean);
psnr_Yl=psnr(Y_low,im_clean);


%% data presentation%%
dictimgh = showdict(D_high,[1 1]*params.blocksize,round(sqrt(params.dictsize)),round(sqrt(params.dictsize)),'lines','highcontrast');
dictimgl = showdict(D_low,[1 1]*params.blocksize,round(sqrt(params.dictsize)),round(sqrt(params.dictsize)),'lines','highcontrast');
figure(); imagesc(imresize(dictimgh,2,'nearest'));colormap('gray');
title('Trained dictionary high');
figure(); imagesc(imresize(dictimgl,2,'nearest'));colormap('gray');
title('Trained dictionary low');

figure();  imshow(Y_high,[]);title(sprintf('Recovered high,  PSNR = %.2fdB', psnr_Yh));
figure();imshow(im_1,[]);title(sprintf('Original,  PSNR = %.2fdB', psnr_Yl));
figure();imshow(im_clean,[]);title(sprintf('Clean Image'));

 

