% function [Y,flag,residual,iter]=iRadon(res1,res2,ErrTol,MaxIts,verbose)
%
% 2-D inverse discrete Radon transform.
% The inverse transform is computed using the conjugate gradient method.
%
%  Input parameters:
%    res1,res22   Discrete Radon sectors as returned from the function Radon.
%    ErrTol       Optional parameter of error tolerance used by the conjugate
%                 gradient method. If not specified, tolerance of 1.e-2 is used.
%    MaxIts       Maximum number of iterations. Default 10.
%    verbose      Display verbose CG information. 0 will suppress verbose information.
%                 Any non-zero value will display verbose CG information.
%
%  Output arguments:
%    Y            The inverted matrix.
%    flag         Convergence flag. See CG for more information.
%    residual     Residual error at the end of the inversion.
%    iter         The iteration number at which ErrTol was achieved. Relevant only if
%                 flag=0
%
% Yoel Shkolnisky 17/12/02

function [Y,flag,residual,iter]=iRadon(res1,res2,ErrTol,MaxIts,verbose)

if nargin<5
   verbose = 0;
end

if nargin<4
   MaxIts = 10;
end

if nargin<3
   ErrTol = 1.e-2;
end

temp1 = cfftd(res1,1);
temp2 = cfftd(res2,1);

[Y,flag,residual,iter]=ippft(temp1,temp2,ErrTol,MaxIts,verbose);

%Revision record
% 15/1/03	Yoel Shkolnisky		Used cfftd instead of column-wise cfft