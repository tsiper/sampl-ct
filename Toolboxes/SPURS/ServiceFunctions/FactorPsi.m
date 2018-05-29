function [L, U, P, Q, R] = FactorPsi( PSI ,ForceFactor,SavePSI )

global PlotFlag DebugFlag;

SPURS_work_dir = getenv('SPURS_DIR');

FullFilePath = [SPURS_work_dir,'/PreCalculated/LU'];
if ~exist(FullFilePath, 'dir')
    mkdir(FullFilePath);
end

HashIn = DataHash(PSI);
FullFileName = [FullFilePath,'/',HashIn,'.mat'];

if (exist(FullFileName, 'file')) && (ForceFactor == 0)
    load(FullFileName)
    disp(['Required PSI was already Factored. Loading PSI factors from ',HashIn,'.mat']);
else
    if (exist(FullFileName, 'file'))
        delete fullFileName;
    end
    [L, U, P, Q, R] = lu (PSI) ;    
    
    if SavePSI == 1
        disp(['Saving PSI factors to file ',HashIn,'.mat']);
        save(FullFileName,'L', 'U', 'P', 'Q', 'R','-v7.3');
    end
end

nnzL = nnz(L);
numelL = numel(L);
fprintf ('L is %d by %d:\n', size(L,1), size(L,2)) ;
fprintf ('L has %d non zeros (out of %d) = %f percent\n', nnzL, numelL,100*nnzL/numelL) ;

if DebugFlag
    plot_title = (['matrix L has density = ',num2str(nnzL/numelL),' with ',num2str(nnzL),' non zeros out of ',num2str(numelL)]);
    h = figure('Name','L matrix','NumberTitle','off');
    spy(L);
    title(plot_title);
end
    
nnzU = nnz(U);
numelU = numel(U);
fprintf ('U is %d by %d:\n', size(U,1), size(U,2)) ;
fprintf ('U has %d non zeros (out of %d) = %f percent\n', nnzU, numelU,100*nnzU/numelU) ;

if DebugFlag
    plot_title = (['matrix U has density = ',num2str(nnzU/numelU),' with ',num2str(nnzU),' non zeros out of ',num2str(numelU)]);
    h = figure('Name','U matrix','NumberTitle','off');
    spy(U);
    title(plot_title);
end