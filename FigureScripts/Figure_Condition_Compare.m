%% Comparing the condition numbers
%   Between the Pseudo-Polar system to that of a regular parallal beam system

Initialize();

% To save or not to save
SaveFlag = 0;

% The resolution grid we are testing on
N_vec = [12:2:32,36:4:64];


%% Running the loop for all the N sizes
cond_rd  = zeros(1,length(N_vec));
cond_rd0 = zeros(1,length(N_vec));
cond_pp  = zeros(1,length(N_vec));
cond_pp0 = zeros(1,length(N_vec));

for i=1:length(N_vec)
    sayhead('Starting Iteration %d, for N=%d',i,N_vec(i));
    % The PB system
    say('Building Parallel-Beam Matrix');
    N = N_vec(i);
    M = 2*N+2;
    N_rd = 2*N+1;
    theta = (0:1/M:1-1/M)*180;
    Hrd = paralleltomo(N,theta,N_rd,N_rd)/(N+1);
    HrdH = full(Hrd'*Hrd);

    % The PP system
    say('Building Pseudo-Polar Matrix');
    Hpp = BuildPPMtx(N);
    HppH = real(Hpp'*Hpp);

    % Preconditioning the PP matrix
    say('Building Preconditioners');
    M_diag = abs(-N:N);
    M_diag(M_diag==0) = 1/N^2;
    M_pp0 = spdiag(repmat(M_diag,[1,2*N+2]));

    % Preconditioning of the PB matrix - The radon transform matrix
    F      = kron(eye(2*N+2),dftmtx(2*N+1))/sqrt(2*N+1); 
    % The filtering matrix
    M_rd0  = real(F'*(spdiag(repmat(fftshift(M_diag),[1,2*N+2])))*F);
    
    % Calculating the Condition Numbers
    say('Calculating Condition Numbers');
    cond_rd(i)  = cond(HrdH);
    cond_rd0(i) = cond(Hrd'*(M_rd0*Hrd));
    cond_pp(i)  = cond(HppH);
    cond_pp0(i) = cond(Hpp'*(M_pp0*Hpp));

end


%% Plotting the condition number graph

% The plot main paramters
SizeX = 1440; SizeY = 960;
TextParams = {'FontSize',27}; 
LineSpec   = {'-.','-','--','-','-','--'};
MarkerSpec = {'+','o','*','d','none','s'};
ColorSpec  = {[.8 0 0],[0 0.5 0],[0 0 .8],[.7 0 .7],[0 0 0],[0 0.5 0.5]};
YLimVec    = [5e-1,1e6]; XLimVec = [10 68];

% The actual plot
figure('Position',[50 50 SizeX,SizeY],'Name','Condition_Number_Compare');
p = semilogy(N_vec,cond_rd,N_vec,cond_rd0,N_vec,cond_pp,N_vec,cond_pp0,N_vec,ones(1,length(N_vec)));
legend('Parallel-Beam','Preconditioned PB','Pseudo-Polar','Preconditioned PP','$\kappa=1$','Location','northwest');

% Styling up the plot
xlim(XLimVec); ylim(YLimVec); set(gca,TextParams{:}); grid on;
xlabel('Number of Detectors - N'); ylabel('Condition Number $\kappa(\cdot)$');
for i=1:length(p)
    p(i).LineStyle = LineSpec{i};     p(i).Marker    = MarkerSpec{i};
    p(i).Color     = ColorSpec{i};    p(i).LineWidth = 2;
    p(i).MarkerSize = 10;
end

% Saving the plot :-)
SaveFigure('Condition_Number_Compare',SizeX,SizeY);