%% Comparing the condition numbers
%   Between the Pseudo-Polar system to that of a regular parallal beam system

sayhead('Outputting figure of Condition Numbers');
Initialize();

N = 8;
N_rd = 2*N+1;
Eps  = 1e-10;

%% Loading a dummy phantom
x = LoadPhantom(N,'zubal');

%% The PB system
say('Building Parallel-Beam Matrix');
M = 2*N+2;
theta = (0:1/M:1-1/M)*180;
Hrd = paralleltomo(N,theta,N_rd,N_rd)/(N+1);
% HrdF = kron(eye(2*N+2),dftmtx(2*N+1))*Hrd;

y_rd = reshape(Hrd*x(:),[N_rd,M]);

HrdH = full(Hrd'*Hrd);

%% The PP system
say('Building Pseudo-Polar Matrix');
Hpp = BuildPPMtx(N);
HppH = real(Hpp'*Hpp);

%% Calculating Condition numbers
say('Calculating Parallel-Beam condition number');
cond_rd = cond(HrdH);
say('Calculating Pseudo-Polar condition number');
cond_pp = cond(HppH);

%% Checking for the ideal preconditioners
say('Finding the best preconditioners');
% M_rd = inv(Hrd*Hrd');
% F = kron(dftmtx(2*N+1),eye(2*N+2))/sqrt(2*N+1);
F  = kron(eye(2*N+2),dftmtx(2*N+1))/sqrt(2*N+1);
F2 = kron(eye(2*N+2),circshift(dftmtx(2*N+1),[N 0]))/sqrt(2*N+1);

% M_pp = inv(Hpp*Hpp');
% M_rd = pinv(full(Hrd'))*pinv(full(Hrd));
M_rd = pinv(Hrd*Hrd'+Eps*eye(size(Hrd,1)));
% M_pp = pinv(Hpp')*pinv(Hpp);
M_pp = pinv(Hpp*Hpp'+Eps*eye(size(Hpp,1)));
% M_pp = 


M_rd0 = spdiag(diag(M_rd));
M_pp0 = spdiag(diag(M_pp));
M_pp1 = spdiag(diag(M_pp));
for i=1:10*N 
    M_pp1 = M_pp1 + sparse(diag(diag(M_pp,i),i)+diag(diag(M_pp,-i),-i));
end

% % Diagonlizing the preconditioners

% F = kron(dftmtx(2*N+2),dftmtx(2*N+1));
% M_rd_fft = F*M_rd;
% % M_pp_fft = F*M_rd;

%% Finding condition numbers with the preconditioners
cond_pp_pc0 = cond(Hpp'*M_pp0*Hpp);
cond_pp_pc_sqrt = cond(Hpp'*(sqrt(real(M_pp0)))*Hpp);
cond_pp_pc1 = cond(Hpp'*M_pp1*Hpp);
codd_pp_pc_real = cond(Hpp'*sqrt(real(M_pp0))*Hpp);
cond_rd_pc = cond(Hrd'*sqrt(M_rd0)*Hrd);

%% Finding the best PP preconditioner according to khatri-rao
HppHppT = real(Hpp*Hpp');
[U, D,V]  = svd(HppHppT);
VdotV = KhatriRao(V,U);
Mopt = pinv(VdotV+Eps*eye(size(VdotV)))*vec(D);
%%
figure;
plot((Mopt));

%% Reconstructing the phantom
x_rd = vec2im(HrdH*x(:));
x_pp = vec2im(HppH*x(:));

x_pp_op = real(App_T(App(x)))*(N+1)^2;

x_rd_pc = real(vec2im(Hrd'*((M_rd*(Hrd*x(:))))));
x_pp_pc = real(vec2im(Hpp'*(M_pp*(Hpp*x(:)))));

% Reconstructing only with a diagonal preconditioner
x_rd0 = real(vec2im(Hrd'*(M_rd0*(Hrd*x(:)))));
x_pp0 = real(vec2im(Hpp'*((M_pp0.^(1/2))*(Hpp*x(:)))));

x_fbp = vec2im(fbp(Hrd,Hrd*(x(:)),theta));
%% Finding the best preconditioner for PP and PB
pp_mat_op_diff = norm(x_pp_op-x_pp,'fro');

%% Plotting results
sayhead('Plotting Results');
say('Plotting the comparison between polar');
figure;
subplot(241); imagesc(x);       title('Ground Truth');
subplot(242); imagesc(x_rd);    title('Parallel Beam'); % colorbar('southoutside');
subplot(243); imagesc(x_rd_pc); title('PB Preconditioner'); %colorbar('southoutside');
subplot(244); imagesc(x_rd0);   title('PB Diagonal'); %colorbar('southoutside');
subplot(245); imagesc(x_pp);    title('Pseudo Polar');  %colorbar('southoutside');
subplot(246); imagesc(x_pp_pc); title('PP Preconditioner');  %colorbar('southoutside');
subplot(247); imagesc(x_pp0);   title('PP Diagonal');  %colorbar('southoutside');
subplot(248); imagesc(x_fbp);   title('FBP');
%%
say('Plotting the comparison between matrices');
figure;
subplot(121);imagesc(abs(M_pp));
subplot(122);imagesc(abs(M_rd));

%% The ammoung of energy on the diagonal
DiagRatioRD = norm(diag(M_rd))/norm(M_rd-diag(diag(M_rd)),'fro');
DiagRatioPP = norm(diag(M_pp))/norm(M_pp-diag(diag(M_pp)),'fro');
DiagRatioPP = sqrt(diag(M_pp)'*diag(M_pp)) / sqrt(sum(vec(abs(M_pp).^2)));
DiagRatioRD = sqrt(diag(M_rd)'*diag(M_rd) / (sum(vec(M_rd.^2))));


%% Plotting the diagonal
figure;
plot(diag(M_rd0)/max(vec(diag(M_rd0))));
hold on;
plot(abs(diag(M_pp0)/max(vec(diag(M_pp0)))));


%% Plotting for paper the diagonal, matrix and others

Width = 1280;
Height = 600;
gap = 0.0;
FontSize = 26;
LineWidth = 2;
mh = [0.02 0.02];
mv = [0.15 0.08];

figure('Position',[50 50 Width, Height]);


subtightplot(1,2,1,gap,mv,mh);imagesc((abs(real(M_pp))));colormap default;
xlabel('$\mathbf{M}_{1}$ - The Ideal PP Preconditioner');
set(gca,'FontSize',FontSize);
c = colorbar;
c.Ticks = [];
set(gca,'YTick',[],'XTick',[]);

MppVec = sqrt(abs(diag(M_pp)));
y_vec = linspace(0,4,4*(2*N+1));
subtightplot(1,2,2,gap,mv,mh);plot(y_vec,MppVec(2*(2*N+1)+1:6*(2*N+1)),'LineWidth',LineWidth);
xlabel('Diagonal of $\mathbf{M}_{1}$ (first 4 repetitions)');
set(gca,'YTick',[],'XTick',[]);

set(gca,'FontSize',FontSize);

SaveFigure('M_pp_Preconditioner',Width,Height);