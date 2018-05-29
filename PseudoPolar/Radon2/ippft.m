% function [Y,flag,residual,iter] = ippft(pp1,pp2,ErrTol,MaxIts,verbose);
%
% Inverse pseudo-polar Fourier transform.
% The inverse transform is computed using conjugate gradient method.
%
% Input arguments:
% pp1,pp2      Pseudo-polar sectors as returned from the function ppft.
% ErrTol       Optional parameter of error tolerance used by the conjugate
%              gradient method. If not specified, tolerance of 1.e-2 is used.
% MaxIts       Maximum number of iterations. Default 10.
% verbose      Display verbose CG information. 0 will suppress verbose information.
%              Any non-zero value will display verbose CG information.
%
% Output arguments:
% Y            The inverted matrix.
% flag         Convergence flag. See CG for more information.
% residual     Residual error at the end of the inversion.
% iter         The iteration number at which ErrTol was achieved. Relevant only if
%              flag=0
%
% Yoel Shkolnisky 18/12/02

function [Y,flag,residual,iter] = ippft(pp1,pp2,ErrTol,MaxIts,verbose);

if nargin<5
   verbose = 0;
end

if nargin<4
   MaxIts = 10;
end

if nargin<3
   ErrTol = 1.e-2;
end

temp = precondAdjPPFT(pp1,pp2);
[Y,flag,residual,iter] = CG('PtP',temp,{},ErrTol,MaxIts,zeros(size(temp)),verbose);
if flag
   warning (sprintf('Inversion did not converge. Residual error %-2.5e',residual));
end
