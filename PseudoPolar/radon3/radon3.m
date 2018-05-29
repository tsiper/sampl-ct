% function rr=Radon3(im)
%
% Fast algorithm for computing the 3-D discrete Radon transform.
% The computation requires O(n^3logn) operations for a nxnxn image.
%
% The algorithm is based on the Fourier slice theorem for the 30D discrete
% Radon transform. 
% The function computes the 3-D pseudo-polar Fourier transform of im and
% then perform 1-D inverse FFT along the parameters k (the pseudo-radius
% parameter).
%
% Due to small round-off errors, the output may contain small imaginary
% components. If the input image is real, the function truncates any
% imaginary components from the output array (since the discrete Radon
% transform of a real image is real)
%
% Input:
%   im - 3D array of size nxnxn.
%        first  index - x direction
%        second index - y direction
%        third  index - z direction
%     In each of the directions, the range 1...n is treated as -n/2...n/2-1.
%     The input image should be of size nxnxn with even n.
%
% Output:
%   xr - 4D array 3x(n+1)x(n+1)x(3n+1) containing the 3-D discrete Radon transform of the input image.
%     The output array has the following format:
%        first  index  - 1 for x-planes, 2 for y-planes, 3 for z-planes
%        second index  - plane intercept
%        third  index  - slope p
%        fourth index  - slope q
%
% Yoel Shkolnisky 03/02/03

function rr=radon3(im)

% verify that the input is a 3D image of size nxnxn
verifyImage(im);

pp = ppft3(im);
rr = icfftd(pp,2);

if isreal(im)
    rr = real(rr);
end
