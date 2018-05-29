%% Intitializing the workspace
Initialize();
DebugFlag = 0;

%% Run Parameters
N            = 64;
DecFactor    = 2; % We have 4*N/DecFactor projection angles
cg_max_iters = 30;
cg_err       = 1e-20;

%% Loading the phantom
x   = LoadPhantom(N,'zubal');

% Designing the decimation factor, to determine how many viewing angles we have
D = DecOperator(DecFactor,'uniform');

% Bulding all the combos and overwriting them accordingly
[ HtH ] = BuildAtA( N ,DecFactor,'uniform' );

%% 
% Generating regular PP measurements
y_pp = Rpp(x);
% Upscaling to work with the fast HtH
y_pp_pad = PP_Upscale(y_pp);

b2 = Hpp_T(y_pp_pad);
x_hat2 = conjgrad(HtH,b2,zeros(size(b2)),cg_max_iters,cg_err);

CompareImages(HtH(b2),Hpp_T(D(Hpp(b2))));
CompareImages(x_hat2,x);
%% Comparing the operators

y_pp_big  = Hpp(x);
b         = Hpp_T(D(y_pp_big));

x_hat = conjgrad(HtH,b,zeros(size(b)),cg_max_iters,cg_err);

%% Plotting
CompareImages(x,x_hat,'Comparing Phantoms');
CompareImages(Hpp(x_hat),y_pp_big,'Comparing Sinograms');

%% Calculating speed gain
iters = 5;
tic;
for i=1:iters
    x1 = HtH(x);
end
t_conv = toc;
tic;
for i=1:iters
    x1 = Hpp_T(D(Hpp(x)));
end
t_old = toc;

time_improvement = t_old / t_conv;