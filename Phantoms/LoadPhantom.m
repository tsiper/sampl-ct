function [ I ] = LoadPhantom( PhantomSize,PhantomType )
%MYPHANTOM Loads up a pre-selected phantom at a specific size

if nargin < 2
    PhantomType = 'shepp';
end

RelativePath = './Phantoms';

N = PhantomSize;

switch PhantomType
    case 'shepp'
        I = phantom(N);
        
    case 'zubal'
        I =reshape(fread(fopen([ RelativePath, '/det_head_u2med.dat' ])),256,256,128);
        SliceNum = 82;
        I = rot90(I(:,:,SliceNum),3);
        I = imresize(I,[N N]);
        I(I<0) = 0;
        I = I./ max(I(:));
        %         I = (I-min(I(:)))./ (max(I(:))-min(I(:)));
        
    case 'brain'
        I = imread([ RelativePath, '/brain.tif' ]);
        I = I(:,:,2);
        I = imresize(double(I),[N N]);
        I = I./max(I(:));
        I(I<0) = 0;
        
        
    case 'thorax'
        I = imread([ RelativePath, '/thorax.tif' ]);
        I = I(21:end-20,21:end-20,2);
        I = imresize(I,[N N]);
        I = double(I)./double(max(I(:)));
        
    case 'sheep'
        Renormalize = @(x) (x-min(x(:)))./range(x(:));
        I=dicomread([ RelativePath, '/sheepscanhigh00003.dcm' ]);
        I=double(I);
        I = Renormalize(I);
        I=imadjust(I,[0.15,0.37]);
        I=imcrop(I,[100 150 256 256]);
        I = imresize(I,[N N]);
        I = Renormalize(I);
        
    case 'circles'
        I = imread([ RelativePath, '/circles.png' ]);
        I = im2double(imresize(I,[N N]));
    case 'more_circles'
        I = imread([ RelativePath, '/more_circles.png' ]);
        I = im2double(imresize(I,[N N]));
    case 'paint'
        I = imread([ RelativePath, '/paint.png' ]);
        I = im2double(imresize(I,[N N]));
    case 'spectral_shepp'
        I = imread([ RelativePath, '/spectral_shepp.png' ]);
        I = im2double(imresize(I,[N N]));
    case 'spectral_zubal'
        I = imread([ RelativePath, '/spectral_zubal.PNG' ]);
        I = im2double(imresize(I,[N N]));
    case 'RLcirc'
        I1 = imread([ RelativePath, '/rightCirc.PNG' ]);
        I1 = im2double(imresize(I1,[N N]));
        I2 = imread([ RelativePath, '/leftCirc.PNG' ]);
        I2 = im2double(imresize(I2,[N N]));
        I = cat(3,I1(:,:,1),I2(:,:,1));
    otherwise
        error(['Please select one of the following phantoms:\n'...
            '''shepp'', ''brain'', ''thorax'', and ''zubal''']);
end

end
