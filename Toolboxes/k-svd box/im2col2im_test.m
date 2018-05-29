
B=imread('images/barbara.png');
B_patch=B(1:16,1:16);
B_vec=B_patch(:);
% C=im2col(B,[1 8]);
% A=col2im(C,[1 8],[size(B,1) size(B,2)],'distinct');
figure(); imshow(B_patch);
title('our first patch');


%% dictionary parameters %%


params.data = B;
params.Tdata = 3;
params.blocksize = 4;
params.dictsize = 64;
params.maxval = 255;
params.trainnum = 40000;
params.iternum = 30;
params.memusage = 'high';
params.sigma=5;

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
params.data = sampgrid(data,blocksize,ids{:});

% remove dc in blocks to conserve memory %
% blocksize = 2000; %%%%????
% for i = 1:blocksize:size(params.data,2)
%   blockids = i : min(i+blocksize-1,size(params.data,2));
%   params.data(:,blockids) = remove_dc(params.data(:,blockids),'columns');
% end

D = ksvd(params,'');
%% dispaly dictionary %%
dictimg = showdict(D,[1 1]*params.blocksize,round(sqrt(params.dictsize)),round(sqrt(params.dictsize)),'lines','highcontrast');
figure(); imshow(imresize(dictimg,2,'nearest'));
title('Trained dictionary');