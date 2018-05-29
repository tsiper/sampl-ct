% function [Y,flag,residual,iter] = iRadon3(res,ErrTol,MaxIts,verbose);
%
% Inverse 3-D discrete Radon transform.
% The inverse transform is computed using conjugate gradient method.
%
% Input arguments:
% res          Discrete 3-D Radon sectors as returned from the function Radon3. res is a 
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
% Yoel Shkolnisky 01/03/03

function [Y,flag,residual,iter] = iradon3(res,ErrTol,MaxIts,verbose)

if nargin<4
   verbose = 0;
end

if nargin<3
   MaxIts = 10;
end

if nargin<2
   ErrTol = 1.e-2;
end

temp = cfftd(res,2);
[Y,flag,residual,iter]=ippft3(temp,ErrTol,MaxIts,verbose);