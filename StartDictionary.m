%% Initialize
Initialize;
global DebugFlag PlotFlag;
DebugFlag = 1;
PlotFlag  = 1;
Renormalize = @(x) (x-min(x(:)))./range(x(:));

%% Simulation Paramters
DecFactor = 4;
DecMethod = 'uniform';
N = 256;
Sigma = 0.005;
L = 50;
MinErrorDiff=0.005;
t = 1;
lambda = 0.1;

%% Load and corrupt image + Radon Transform
x0 = LoadPhantom(N,'sheep');

y0 = invF1(App(x0));

y = y0 + Sigma*range(y0(:))*randn(size(y0));
y = DecFactor*D(F1(y),DecFactor,DecMethod);

%% Load Dictionaries
DLow = load('d:/DLow211116.mat');
DHigh = load('d:/DHigh211116.mat');
% DLow = load('d:/DLow211116PS25.mat');
% DHigh = load('d:/DHigh211116PS25.mat');
DLow = DLow.D_low;
DHigh = DHigh.D_high;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  START Reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = zeros(size(x0));
u = zeros(size(x0));

iters = 20;
ii = iters;
NumIter=size(z,2);
PARAMS.blocksize=sqrt(size(DLow,1));
PARAMS.maxatoms=PARAMS.blocksize-1;
PARAMS.OverlapSize=PARAMS.blocksize-3;
PARAMS.MinErrDiff=MinErrorDiff;

G=DLow'*DLow;
figure;
cpb = ConsoleProgressBar;
cpb.setText('Dictionary CT Reconstruction');
cpb.start();
while ii > 0
    z = z - t/L*App_T(M(App(u)-y));
    u = FGP(z,2*lambda/L,50);
    [SuperRes,a]=BuildSuperRes_ct(u,DLow,DHigh,PARAMS);
%     [SuperRes,a]=BuildSuperRes_ct(z,DLow,DHigh,PARAMS);
%     u = SuperRes;
    t = (1+sqrt(1+4*t^2))/2;
    ii = ii -1;
    %% TODO Add FISTA acceleration step
    
    if DebugFlag 
        subplot(131);imagesc(z);
        subplot(132);imagesc(u);
        subplot(133);imagesc(SuperRes);
        drawnow;
    end
    cpb.setValue((iters-ii)/iters);
end
cpb.stop();
% SuperRes = u;
% SuperRes=FGP(SuperRes,2*lambda/L,10);
% for i=1:size(a,2)
%     s(i)=length(nonzeros(a(:,i)));
% end
Mean_coeff=mean(sum(a~=0));



%% Plotting
if PlotFlag
    Ah=1;
    Al=1;
    u=RenormalizeImage(u);
    x0=RenormalizeImage(x0);
    xBad=Renormalize(App_T(M(y)));
    SuperRes=Renormalize(SuperRes);
    Ah=fminsearch(@(Ah)-ssim(u,x0*Ah),1);
    Al=fminsearch(@(Al)-ssim(u,x0*Al),1);
    SsimH=ssim(u,Ah*x0);
    SsimL=ssim(xBad,Al*x0);
    CompareImages(xBad,x0);
    CompareImages(SuperRes,x0);
    CompareImages(u,x0);
end
