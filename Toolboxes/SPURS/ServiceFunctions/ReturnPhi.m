function [ PHI ] = ReturnPhi( InputGridCoordinates, OutputGridCoordinates, KernelFunction, ForceGenrateNew )

PHI = NaN;
PHIsettings.KernelFunction = KernelFunction;
PHIsettings.InputGridCoordinates = InputGridCoordinates;
PHIsettings.OutputGridCoordinates = OutputGridCoordinates;

SPURS_work_dir = getenv('SPURS_DIR');
HashIn = DataHash(PHIsettings);
FullFilePath = [SPURS_work_dir,'/PreCalculated/PHI'];
if ~exist(FullFilePath, 'dir')
    mkdir(FullFilePath);
end
FullFileName = [FullFilePath,'/',HashIn,'.mat'];


if (exist(FullFileName, 'file')) && (ForceGenrateNew == 0)
    load(FullFileName)
    disp(['Required PHI was already computed. Loading PHI from ',HashIn,'.mat']);
else
    if (exist(FullFileName, 'file'))
        delete fullFileName;
    end
    
    switch(KernelFunction.type)
        case 'Bspline'
            disp(['Building PHI for a ',KernelFunction.type,' kernel']);
            r = InputGridCoordinates(:,1) + 1j.*InputGridCoordinates(:,2);
            t = OutputGridCoordinates(:,1) + 1j.*OutputGridCoordinates(:,2);
            PHI=BsplinePHI(r, t,KernelFunction.degree);
            disp(['Saving PHI to file ',HashIn,'.mat']);
            save(FullFileName,'PHI','-v7.3');
        otherwise
            disp('ERROR - kernel ',KernelFunction.type,' not supported');
            return;
    end
    
    
end
nnzPHI = nnz(PHI);
numelPHI = numel(PHI);
fprintf ('PHI is %d by %d:\n', size(PHI,1), size(PHI,2)) ;
fprintf ('PHI has %d non zeros (out of %d) = %f percent\n', nnzPHI, numelPHI,100*nnzPHI/numelPHI) ;
end

