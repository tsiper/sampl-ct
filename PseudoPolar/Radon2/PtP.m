% function Y = PtP(X)
%
% Gram Operator of the pseudo-polar Fourier transform.
% Performs adjP(D(P)) where
%    P     The pseudo-polar Fourier transform
%    D     Preconditioner
%    adjP  Adjoint pseudo-polar Fourier transform
%
%  Input parameters:
%    X      n*n matrix (x,y)
%  Outputs parameters:
%    Y      n*n matrix (x,y)
%
% Yoel Shkolnisky 17/12/02

function Y = PtP(X)
[pp1,pp2] = OptimizedPPFT(X);
Y = precondAdjPPFT(pp1,pp2);
