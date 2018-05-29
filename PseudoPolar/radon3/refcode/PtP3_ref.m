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
% Yoel Shkolnisky 28/2/03
%
% Revisions:
% Yoel Shkolnisky 19/05/2013    Renamed from PtP3 to PtP3_ref.

function Y = PtP3_ref(X)
pp = ppft3_ref(X);
Y = precondadjppft3_ref(pp);
