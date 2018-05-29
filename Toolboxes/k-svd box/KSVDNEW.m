
%input data%
im_in=double(imread('images/barbara.png'));
im_in=im_in(:,:,1);
params.data = im_in;
params.Tdata = 3;
params.blocksize = 16;
params.dictsize = 256;
params.maxval = 255;
params.trainnum = 100000;
params.iternum = 15;
params.memusage = 'high';
%params.sigma=5;%
D=ksvd(params,'');
Imtcol=im2colstep(im_in,[params.blocksize, params.blocksize]);
SparseMat=omp(D'*Imtcol,D'*D,params.Tdata);
Y= col2imstep(Y,[size(im_in,1) size(im_in,1)],[params.blocksize params.blocksize]);
dictimg = showdict(D,[1 1]*params.blocksize,round(sqrt(params.dictsize)),round(sqrt(params.dictsize)),'lines','highcontrast');
figure(); imshow(imresize(dictimg,2,'nearest'));
title('Trained dictionary');
 figure; imshow(Y);


