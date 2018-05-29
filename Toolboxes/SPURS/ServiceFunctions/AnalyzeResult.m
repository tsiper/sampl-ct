function [ measures ] = AnalyzeResult(ResultImages,RefImg,NameString)

sqrtN = size(ResultImages,1);
ResultImages = real(ResultImages);
RefImg = real(RefImg);
for ii=1:size(ResultImages,3)
    ResultImage = ResultImages(:,:,ii);
    
    UseCircMaskOnResult = 1;
    if (UseCircMaskOnResult)
        Nd2 = sqrtN/2;
        [xx, yy] = meshgrid(-Nd2:Nd2-1,-Nd2:Nd2-1);
        circlemask = sqrt(xx.^2+yy.^2)<=Nd2;
        ResultImage = circlemask.*ResultImage;
        RefImg = circlemask.*RefImg;
        gain = (ResultImage(:)'*RefImg(:))/(ResultImage(:)'*ResultImage(:));
        ResultImage = ResultImage.*gain;
    end
    
    MeasuresStrings = {'MSE','MAXE','RMSE','SNR','PSNR','MSSIM'};
    CalcDB = [0,0,0,0,0,0];
    MeasuresSuffix = {'','','',' dB',' dB',''};
    
    measures.MSE(ii) = mean((ResultImage(:)-RefImg(:)).^2);
    measures.MAXE(ii) = max(abs(ResultImage(:)-RefImg(:)));
    measures.RMSE(ii) = sqrt(measures.MSE(ii));
    [measures.PSNR(ii),measures.SNR(ii)] = psnr(ResultImage,RefImg);
    [measures.MSSIM(ii),~] = ssim(ResultImage,RefImg);
    
end

disp(['Image quality measures for ',NameString]);
for ii=1:size(MeasuresStrings,2)
    MeasureString = cell2mat(MeasuresStrings(:,ii));
    MeasureSuffix = cell2mat(MeasuresSuffix(:,ii));
    if ii<=3
        [BestV BestI] = min(getfield(measures,MeasureString));
    else
        [BestV BestI] = max(getfield(measures,MeasureString));
    end
    if(CalcDB(ii))
        BestV = db(BestV);
    end
    disp(['Best ',MeasureString,' is ',num2str(BestV),MeasureSuffix,', for iteration ',num2str(BestI)]);
    %     FigHandle = figure(10+ii);
    FigHandle = figure('Name',MeasureString,'NumberTitle','off');
    
    subplot(1,3,1);
    imagesc(ResultImages(:,:,BestI));axis square; colormap('gray'); axis off;%colorbar
    title(['Image with best ',MeasureString,'=',num2str(BestV),MeasureSuffix,', @ it.',num2str(BestI)]);
    
    subplot(1,3,2);
    plot(1:size(ResultImages,3),getfield(measures,MeasureString));
    xlabel('Iteration #');
    ylabel([MeasureString,' ',MeasureSuffix]);
    title([NameString,' ',MeasureString,' measure']);
    set(FigHandle,'Position',[50 50 1500 400],'PaperPositionMode','auto');
    
    subplot(1,3,3);
    wanted_line = round(113*sqrtN/256);
    plot((-sqrtN/2:sqrtN/2-1),ResultImages(wanted_line,:,BestI),'-','LineWidth',2,'Color',[0,0,0]); hold on;
    plot((-sqrtN/2:sqrtN/2-1),RefImg(wanted_line,:),'--','LineWidth',1,'Color',[0.4660,0.6740,0.1880]); hold off;
    grid on;ylim([-0.1 1.3])
    title(['Line ',num2str(wanted_line),' of image with best ',MeasureString,'.']);
    
    
end

end

