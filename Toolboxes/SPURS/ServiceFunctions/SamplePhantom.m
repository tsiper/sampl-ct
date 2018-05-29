function [RefImg, PhantomSamples] = SamplePhantom(PhantomString,TrajrctoryString,sqrtN,M,GridCoordinates,ForceCalc)

N = sqrtN^2;

FileName = [PhantomString,'_',TrajrctoryString,'_N',num2str(sqrtN^2),'_M',num2str(M)];
SPURS_work_dir = getenv('SPURS_DIR');
FullFilePath = [SPURS_work_dir,'/PreCalculated/PhantomData'];

if ~exist(FullFilePath, 'dir')
    mkdir(FullFilePath);
end

FullFileName = [FullFilePath,'/',FileName,'.mat'];

if (exist(FullFileName, 'file')) && (ForceCalc == 0)
    load(FullFileName)
    disp(['Required analytical Phantom samples were already calaulated. Loading Phantom data from ',FileName,'.mat']);
else
    if (ForceCalc == 1)
        disp(['ForceCalc flag = 1, therefore starting calculation of analytical Phantom sample values']);        
    else
        disp([FileName,'.mat was not found. Starting calculation of analytical Phantom sample values']);
    end
    switch PhantomString
        case 'Phantom1'
            InputSamplingPointCoor = [GridCoordinates(:,1) GridCoordinates(:,2)].';
            sqrmat = makephantom1(sqrtN,1);
            [RefImg, simfid] = squares_fast(sqrtN,sqrmat,InputSamplingPointCoor,1,0);
            simfid = simfid./N;
            RefImg = rot90(RefImg.',2);
            PhantomSamples = [GridCoordinates(1:M,1) GridCoordinates(1:M,2) real(simfid.') imag(simfid.')];
            
        case 'AnalyticalSL'
            InputSamplingPointCoor = [GridCoordinates(1:M,2) GridCoordinates(1:M,1)].';
            [RefImg, simfid] = MakeNewPhantom(sqrtN,InputSamplingPointCoor,'SL');
            PhantomSamples = [GridCoordinates(1:M,1) GridCoordinates(1:M,2) real(simfid.') imag(simfid.')];

            
        case 'Brain'
            InputSamplingPointCoor = [GridCoordinates(1:M,2) GridCoordinates(1:M,1)].';
            [RefImg, simfid] = MakeNewPhantom(sqrtN,InputSamplingPointCoor,'Brain');     
            PhantomSamples = [GridCoordinates(1:M,1) GridCoordinates(1:M,2) real(simfid.') imag(simfid.')];

    end
    if (exist(FullFileName, 'file'))
        Old=load(FullFileName,'PhantomSamples');
        if sum(abs(Old.PhantomSamples(:) - PhantomSamples(:))) == 0
            disp(['The new generated samples are the same as the ones in ',FileName,'.mat']);
        else
            disp(['The new generated samples are different from the ones in ',FileName,'.mat']);
        end
        disp(['Deleting old data file ',FileName,'.mat']);
        oldFolder = cd(FullFilePath);
        delete([FileName,'.mat']);
        cd(oldFolder);
    end
    disp(['Saving analytical Phantom data to file ',FileName,'.mat']);
    save(FullFileName,'PhantomSamples','RefImg','-v7.3');
end
end

