% function [Y,flag,residual,iter] = ippft3(pp,ErrTol,MaxIts,verbose);
%
% Inverse 3-D pseudo-polar Fourier transform.
% The inverse transform is computed using conjugate gradient method.
%
% Input arguments:
% pp           Pseudo-polar samples as returned from the function ppft3. pp is a 
%              4-D array of size 3x(3n+1)x(n+1)x(n+1).
% ErrTol       (Optional) Error tolerance used by the conjugate
%              gradient method. Default 1.e-2.
% MaxIts       Maximum number of iterations. Default 10.
% verbose      Display verbose CG information. 0 will suppress verbose information.
%              Any non-zero value will display verbose CG information.
%
% Output arguments:
% Y            The reconstructed 3-D array.
% flag         Convergence flag. See CG for more information.
% residual     Residual error at the end of the inversion.
% iter         The iteration number at which ErrTol was achieved. Relevant only if
%              flag=0
%
% Yoel Shkolnisky 26/02/03
%
% Revision:
% Yoel Shkolnisky 19/05/2013 OptimizedprecondAdjPPFT3 was renamed to
%      precondAdjPPFT3. Replaced call accordingly.


function [Y,flag,residual,iter] = ippft3(pp,ErrTol,MaxIts,verbose)

if nargin<4
   verbose = 0;
end

if nargin<3
   MaxIts = 10;
end

if nargin<2
   ErrTol = 1.e-2;
end

temp = precondadjppft3(pp);
[Y,flag,residual,iter] = CG('PtP3',temp,{},ErrTol,MaxIts,zeros(size(temp)),verbose);
if flag
   warning ('Inversion did not converge. Residual error %-2.5e',residual);
end
