function [ ] = CompareImages( y,y_ref,Title,varargin )
%COMPARESINOGRAM Compares two sinograms, with all the relevant plots and other
%useful information

[m,n] = size(y);
[M,N] = size(y_ref);

if nargin<3
    Title = '';
end

Eps = 1e-3;
myfft = @(x) db(abs(fftshift(fft2(x,4*size(x,1),4*size(x,2))))+Eps);

% Setting default flags
fftFlag = 1;
rows    = 2;
titleFlag = 1;
FigWidth = 1280; 
FigHeight = 720;

for k=1:length(varargin)
    if strcmpi(varargin{k},'nofft')
        fftFlag = 0;
        rows    = 1;
        FigHeight = 600;
    elseif strcmpi(varargin{k},'nodb')
        myfft = @(x) (abs(fftshift(fft2(x,4*size(x,1),4*size(x,2))))+Eps);
    elseif strcmpi(varargin{k},'notitle')
        titleFlag = 0;
    end
end

if (m~=M)||(n~=N)
    error('Images aren''t the same size');
end

%% Calcaulting Statistics
PSNR = psnr(y,y_ref);
RMSE = rmse(y,y_ref);

% PSNRshift = psnr2(y,y_ref);

% Find Amplitude offset
% AmpOffset = fminsearch(@(A) -psnr(A*y,y_ref),1);

% Find registration offset
% shift = FindShift( y,y_ref );

%% Plotting results
gap = 0.1;
mh = 0.05;
mw = [0.03 0.08];
figure('Name',Title,'Position',[50 -50 FigWidth FigHeight]);
ax(1) = subtightplot(rows,3,1,gap,mh,mw); imagesc(y_ref);colormap gray;
title('Ground Truth');
ax(2) = subtightplot(rows,3,2,gap,mh,mw); imagesc(y);colormap gray;
% title(['PSNR=',num2str(PSNR),' PSNRshift=',num2str(PSNRshift)]);
title(['PSNR=',num2str(PSNR)]);
ax(3) = subtightplot(rows,3,3,gap,mh,mw); imagesc(y-y_ref);colormap gray;
title(['Diff Image RMSE=',num2str(RMSE)]);
if fftFlag
    ax2(1) = subtightplot(rows,3,4,gap,mh,mw); imagesc(myfft(y_ref));colormap gray;
    title('Original fft');
    ax2(2) = subtightplot(rows,3,5,gap,mh,mw); imagesc(myfft(y));colormap gray;
    title('Recon. fft');
    ax2(3) = subtightplot(rows,3,6,gap,mh,mw); imagesc(myfft(y)-myfft(y_ref));colormap gray;
    title('fft diff');
end
if titleFlag
    suptitle(Title);% ,' Offset=[',num2str(shift(:)'),']']);
else
    suptitle(['Amplitude=',num2str(AmpOffset)]);% ,' Offset=[',num2str(shift(:)'),']']);
end
MyColorbar(ax); 
linkaxes(ax);   
if fftFlag
    linkaxes(ax2);
    MyColorbar(ax2);
end


end

