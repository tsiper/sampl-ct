N = 10;
PolarGrid = BuildPolarGrid(N,8);
PolarGridud = flipud(PolarGrid); % Polar grid, built backwards
ang = 270;  
rot_matrix= [cosd(-ang) -1*sind(-ang) ; sind(-ang) cosd(-ang)]; %around z
PolarGridRot = (rot_matrix*PolarGrid');
PolarGridRot = PolarGridRot';
figure;
% for i=1:length(PolarGridRot)
%     scatter(PolarGridRot(i,1),PolarGridRot(i,2));
%     xlim([-N/2-1,N/2+1]);
%     ylim([-N/2-1,N/2+1]);
%     hold on;
%     pause(0.1);
% end

%% 
N = 128;
Mp = 180;
x0 = LoadPhantom(N,'zubal');
N_rd = 2*ceil(norm(size(x0)-floor((size(x0)-1)/2)-1))+3;
theta = (0:1/Mp:1-1/Mp)*180;

y = radon(x0,theta);
figure;subplot(121);imagesc(x0);subplot(122);imagesc(y);
sFactor = N/N_rd;
PolarGrid = BuildPolarGrid(N_rd,Mp)*sFactor;

%%
SPURS_settings.sqrtN = N;
SPURS_settings.KernelFunctionString = 'Bspline';
SPURS_settings.KernelFunctionDegree = 3;
SPURS_settings.ReusePrecalculatedData = 0;
SPURS_settings.Rho = 1e-3;
SPURS_settings.Niterations = 5;
SPURS_settings.UseW = 0;
SPURS_settings.ForceGenrateNewPhi = 0;
SPURS_settings.ForceFactorPsi = 0;
SPURS_settings.SavePSI = 0;
SPURS_settings.OverGridFactor = 1;
SPURS_settings.alpha = 1;
SPURS_settings.CalcOptimalAlpha = 1;
SPURS_settings.FilterInImageSpace = 1;

%% Running SPURS
Samples = -vec(F1(y)/N^2);
[ OutputImages, b_hat] = SPURS(Samples, PolarGrid, SPURS_settings);

%% Running Filtered BackProjection
x_fbp = iradon(y,theta);
x_fbp = x_fbp(2:end-1,2:end-1);
x_fbp2= EqualizeImage(ShiftImage(x_fbp,x0),x0);

%% Shift Images
x_hat = OutputImages(:,:,1);
x_hat = ShiftImage(x_hat,x0);
x_hat2 = EqualizeImage(x_hat,x0);

%%
CompareImages(x_hat,x0,'HAT');
CompareImages(x_hat2,x0,'Equalized HAT2');
CompareImages(x_hat,x_hat2,'Hat1 Vs. Hat2');
CompareImages(x_fbp,x0);