function [ SPURS_settings ] = SPURS_DefaultParams( N )
%SPURS_DEFAULTPARAMS Returns a basic set of SPURS settings

% Setting up the enviromental parameters
setenv('SPURS_DIR', pwd);
setenv('SPURS_RUN_TIME', datestr(now,'ddmmyyyyTHHMMSS'));

% The default settings
SPURS_settings.sqrtN = N;
SPURS_settings.KernelFunctionString = 'Bspline';
SPURS_settings.KernelFunctionDegree = 3;
SPURS_settings.ReusePrecalculatedData = 1;
SPURS_settings.Rho = 1e-3;
SPURS_settings.Niterations = 5;
SPURS_settings.UseW = 0;
SPURS_settings.ForceGenrateNewPhi = 0;
SPURS_settings.ForceFactorPsi = 0;
SPURS_settings.SavePSI = 1;
SPURS_settings.OverGridFactor = 1;
SPURS_settings.alpha = 1;
SPURS_settings.CalcOptimalAlpha = 1;
SPURS_settings.FilterInImageSpace = 1;

end

