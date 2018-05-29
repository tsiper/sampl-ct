% function [r,phi,theta] = pp2spherical(sector,k,l,j,n)
%
% Convert pseudo-polar coordinates to spherical coordinates.
%
% Convert the cartesian coordinates of the pseudo-polar sample 
% pp(sector,k,l,j) into spherical coordinates expressed in terms of
% (r,phi,theta).
% The cartesian coordinates of the sample pp(sector,k,l,j) are given by
%    sector=1: pp(1,k,l,j) is the frequnecy sample at (k,-2lk/n,-2jk/n)
%    sector=2: pp(2,k,l,j) is the frequnecy sample at (-2lk/n,k,-2jk/n)
%    sector=3: pp(3,k,l,j) is the frequnecy sample at (-2lk/n,-2jk/n,k) 
% The function accepts the indices (sector,k,l,j), computes the
% corresponding cartesian sample and converts this cartesian sample into
% spherical sample.
%
% Parameters:
%       sector    Index of the pseudo-polar sector. must be 1, 2, or 3.
%       k,l,j     Index of the pseudo-polar sample.
%       n         Side length of the original image.
%
% See also ppft3
%
% Yoel Shkolnisky 03/03/03

function [r,phi,theta] = pp2spherical(sector,k,l,j,n)

% convert k,l,j into cartesian coordinates
if (sector==1)
    x = k;
    y = -2*l*k/n;
    z = -2*j*k/n;
elseif (sector==2)
    x = -2*l*k/n;
    y = k;
    z = -2*j*k/n;
elseif (sector==3)
    x = -2*l*k/n;
    y = -2*j*k/n;
    z = k;
else
    error('sector must be 1,2, or 3');
end

% convert the cartesian coordinates into spherical coordinates
r     = sqrt(x^2+y^2+z^2);
phi   = acos(z/r);
theta = atan(y/x);

    