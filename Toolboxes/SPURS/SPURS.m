function [ OutputImages, b_hat] = SPURS(b, Kappa_m, SPURS_settings)

global PlotFlag DebugFlag;

b_hat = NaN;
M = length(b);
sqrtN = SPURS_settings.sqrtN;

if SPURS_settings.OverGridFactor ~= ceil(sqrtN*SPURS_settings.OverGridFactor/2)*2/sqrtN
    disp(['Over grid factor was corrected from ',num2str(SPURS_settings.OverGridFactor),' to ',num2str(ceil(sqrtN*SPURS_settings.OverGridFactor/2)*2/sqrtN)]);
    SPURS_settings.OverGridFactor = ceil(sqrtN*SPURS_settings.OverGridFactor/2)*2/sqrtN;
end

OverGridFactor = SPURS_settings.OverGridFactor;

ReconstructionGridCoordinates = ConstructGridCartesian( SPURS_settings.sqrtN, 1, 0, 0, SPURS_settings.OverGridFactor);

KernelFunction.type = SPURS_settings.KernelFunctionString;
KernelFunction.degree = SPURS_settings.KernelFunctionDegree;

[ PHI ] = ReturnPhi( ReconstructionGridCoordinates, Kappa_m, KernelFunction, SPURS_settings.ForceGenrateNewPhi);
[m, n] = size(PHI);
W = speye(M);
[ PSI ] = ReturnPsi( PHI, W, SPURS_settings.Rho, SPURS_settings.UseW );
clear PHI;

[L, U, P, Q, R] = FactorPsi( PSI ,SPURS_settings.ForceFactorPsi,SPURS_settings.SavePSI );

clear PSI;

[ r_aq,ImgFilter ] = ReturnFilters( KernelFunction,OverGridFactor,sqrtN );

disp(['SPURS - Using original samples (b) RMSe (iteration ',num2str(1),')']);
if (SPURS_settings.FilterInImageSpace == 1)
    disp(['SPURS - Filtering is performed in image-space']);
else
    disp(['SPURS - Filtering is performed in k-space']);
end
disp(['SPURS - Starting first iteration']);
b_n = b;

OutputImages = zeros(sqrtN,sqrtN,SPURS_settings.Niterations);
samplesRMSE = NaN;

for ii = 1:SPURS_settings.Niterations
    if(SPURS_settings.UseW)
        b_tag = [W*b_n ; zeros(n,1)];
    else
        b_tag = [b_n ; zeros(n,1)];
    end
    
    c_tag = Q * (U \ (L \ (P * (R \ b_tag)))) ;
    c = c_tag(m+1:end);
    c = reshape(c,sqrt(n),sqrt(n));
    
    if (SPURS_settings.FilterInImageSpace == 1)
        [uniform_k_samples, OutputImage] = IFFTandFilterInImageSpace(c,ImgFilter,sqrtN,OverGridFactor);
    else
        [uniform_k_samples, OutputImage] = FilterInKSpaceandIFFT(c,r_aq,sqrtN,OverGridFactor);
    end
    
    if ((ii==1) || (ii==SPURS_settings.Niterations)) && DebugFlag
        if(ii==1)
            figure('Name','SPURS 1st It.','NumberTitle','off');
        else
            figure('Name',['SPURS It. ',num2str(ii),' '],'NumberTitle','off');
        end
        imagesc(real(OutputImage));axis square; colormap('gray'); axis off;
        title(['SPURS Result after ',num2str(ii),' iterations']);
    end
    OutputImages(:,:,ii) = OutputImage;
    
    if (SPURS_settings.Niterations == 1)
        continue;
    end
    [ b_hat ] = Interpolate_with_sinc(sqrtN, ReconstructionGridCoordinates, uniform_k_samples(:), Kappa_m );
    sample_err = b-b_hat;
    samplesRMSE(ii) = sqrt(mean((abs(sample_err)).^2));
    disp(['SPURS - Resampling error at sampling coordinates: RMSE(b_0-b_',num2str(ii),')=',num2str((samplesRMSE(ii))),' (',num2str(db(samplesRMSE(ii))),' [dB])']);
    
    
    if (SPURS_settings.CalcOptimalAlpha)
        sample_err_tag = [sample_err ; zeros(n,1)];
        sample_err_c_tag = Q * (U \ (L \ (P * (R \ sample_err_tag)))) ;
        sample_err_c = sample_err_c_tag(m+1:end);
        sample_err_c = reshape(sample_err_c,sqrt(n),sqrt(n));
        
        if (SPURS_settings.FilterInImageSpace == 1)
            [sample_err_uniform_k_samples, ~] = IFFTandFilterInImageSpace(sample_err_c,ImgFilter,sqrtN,OverGridFactor);
        else
            [sample_err_uniform_k_samples, ~] = FilterInKSpaceandIFFT(sample_err_c,r_aq,sqrtN,OverGridFactor);
        end
        
        [ sample_err_b_hat ] = Interpolate_with_sinc(sqrtN, ReconstructionGridCoordinates, sample_err_uniform_k_samples(:), Kappa_m );
        alpha = real(sample_err'*sample_err_b_hat/sum((abs(sample_err_b_hat)).^2));
    end
    b_n = b_n + alpha.*sample_err;
end
if (SPURS_settings.Niterations > 1)&&(DebugFlag||PlotFlag)
    FigHandle = figure('Name','SPURS - Error on sampling locations','NumberTitle','off');
    plot(1:SPURS_settings.Niterations,db(samplesRMSE));
    xlabel('Iteration #');
    ylabel('RMSE(b_0-b_n) [dB]');
    title(['SPURS - Error on sampling locations RMSE(b_0-b_n)']);
    set(FigHandle,'Position',[100 100 500 400],'PaperPositionMode','auto');
end