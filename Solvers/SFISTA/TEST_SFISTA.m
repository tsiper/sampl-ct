clc;
clear;
close all;

% This script tests the Smoothing based FISTA algorithm for solving the analysis sparse recovery problem
% -------------------------------------------------------------------------------------------------------
rng('default'); rng(1);         % For reproducibility

%% Parameters
% -------------------------------------------------------------------------
M       = 70; 
N       = 128;
P       = 144;
K       = 15;

Nsig    = 0;                               % Noise STD

% Create experiment
A       = randn(M, N);                  % Gaussianly random sensing matrix

% Normalize the columns of A
A       = NormMat( A, 'l2' );

% Support in the sparse domain
Supp    = randperm(N, K);               % Random support

z       = zeros(N, 1);
z(Supp) = randn(length(Supp), 1);       % Gaussianly random K-sparse vector

% -------------------------------------------------------------------------
% Generate x - dense vector, but sparse in some other basis:
% x = D*z, z is K-sparse 
% At this stage I use only unitary transformations, such that DD^* = I, That is D^* = D^{-1}
% -------------------------------------------------------------------------
TransformType = 'fwht';                   % Type of transformation
switch lower(TransformType)
    case 'wave' % Wavelet transform - based on the Rice Wavelet Toolbox, Version 2.4, Dec 2002
        WaveD     = 3;                                % Wavelet depth hirarchy
        WaveScale = daubcqf(16,'mid');                % Type of wavelet filter
        D         = @(x) mdwt(x, WaveScale, WaveD);   % Synthesis operator
        Dt        = @(x) midwt(x, WaveScale, WaveD);  % Analysis operator
        x  = D(z);
        
        Dnorm = 0.1;
    case 'fwht' % Fast Welch Hadamard transform
        D  = @fwht;                                   % Synthesis operator
        Dt = @ifwht;                                  % Analysis operator
        x  = D(z);
        
        Dnorm = 0.1;
    case 'dct'  % DCT
        D  = @dct;                                    % Synthesis operator
        Dt = @idct;                                   % Analysis operator
        x = D(z);
        
        Dnorm = norm(dctmtx(N), 2);
    case 'qr'   % Tight frame usinq QR factorization
        D     = randn(N, P);
        [Q,R] = qr(D);
        Dt    = Q(:, 1:N);
        D     = Dt';                                  % This results in a Unitary D, so inv(D) = D' (D is l2-normalized)
        x     = D*z;
        
        Dnorm = norm(D, 2);
end

% Generate data
b = A*x + Nsig*randn(M, 1);

% SFISTA parameters
SFISTA_Params.Lambda       = 4e-3;
SFISTA_Params.mu0          = 1e-1;
SFISTA_Params.muf          = 1e-3;
SFISTA_Params.mu           = (1e-3)/SFISTA_Params.Lambda;
SFISTA_Params.gamma        = 3;   
SFISTA_Params.Continuation = 0;
SFISTA_Params.IterMax      = 1000;
SFISTA_Params.L            = norm(A, 2)^2 + (Dnorm^2)/SFISTA_Params.mu;        % If D is unitary then norm(D, 2) = 1, I think

%% Solve the problem
% -------------------------------------------------------------------------
% Solve SFISTA
tic;
[x_sfista, F_stack] = SFISTA( b, A, [], D, Dt, SFISTA_Params );
toc;

%% Plot the result
% -------------------------------------------------------------------------
figure;
loglog(F_stack, '-*'); grid on; title('Function value in log scale');

figure;
subplot(2,1,1);
stem(1:N, x);title(['x vector, transform type: ' TransformType]); hold on; grid on;
stem(1:N, x_sfista, 'm+'); legend('x', 'x_{sfista}');
subplot(2,1,2);
stem(1:N, z);title('z vector'); hold on; grid on;
switch lower(TransformType)
    case 'qr'
        stem(1:N, Dt*x_sfista, 'm+'); legend('z', 'z_{sfista}');
    otherwise
        stem(1:N, Dt(x_sfista), 'm+'); legend('z', 'z_{sfista}');
end



