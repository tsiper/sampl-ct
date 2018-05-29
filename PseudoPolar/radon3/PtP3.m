% function Y = PtP3(X)
%
% Gram Operator of the 3-D pseudo-polar Fourier transform.
% Performs adjP(D(P)) where
%    P     The 3-D pseudo-polar Fourier transform
%    D     Preconditioner
%    adjP  Adjoint 3-D pseudo-polar Fourier transform
%
%  Input parameters:
%    X      nxnxn 3-D array (n even)
%  Outputs parameters:
%    Y      nxnxn 3-D array
%
% Yoel Shkolnisky, July 2007.

function Y = PtP3(X)
pp = ppft3(X);
Y = precondadjppft3(pp);
