function [ PSI ] = ReturnPsi( PHI, W, Rho, UseW )

global DebugFlag;

[m, n] = size(PHI);

if(UseW)
    PSI = [speye(m) W*PHI; (W*PHI).' -(Rho).*speye(n)];
else
    PSI = [speye(m) PHI; PHI' -(Rho).*speye(n)];
end
nnzPSI = nnz(PSI);
numelPSI = numel(PSI);

if DebugFlag
    plot_title = (['matrix PSI LS has density = ',num2str(nnzPSI/numelPSI),' with ',num2str(nnzPSI),' non zeros out of ',num2str(numelPSI)]);
    h = figure('Name','PSI matrix','NumberTitle','off');
    spy(PSI);
    title(plot_title);
end
fprintf ('PSI is %d by %d:\n', size(PSI,1), size(PSI,2)) ;
fprintf ('PSI has %d non zeros (out of %d) = %f percent\n', nnzPSI, numelPSI,100*nnzPSI/numelPSI) ;

end

