% function [res1,res2]=Radon(im);
%
% Fast algorithm for computing the discrete Radon transform.
% The computation requires O(n^2logn) operations.
%
% im    The image whose discrete Radon transform should be computed.
%       Must be real (no imaginary components) and of a dyadic square size (2^k x 2^k).
%
% Returns res1 and res2 (of size 2n+1xn+1) that contain the discrete Radon
% transform of the input image im.
% res1 contains the values which correspond to rays of radius k=-n...n and angles
% theta=arctan(2l/n) l=-n/2...n/2 (from -pi/4 to pi/4). 
% res2 contains the values which correspond to rays of radius k=-n...n and angles
% theta=pi/2-arctan(2l/n) l=-n/2...n/2 (from 3pi/4 to pi/4). 
%
% Due to small round-off errors, the output may contain small imaginary
% components. If the input image is real, the function truncates any
% imaginary components from the output array (since the discrete Radon
% transform of a real image is real)
%
% See thesis' final PDF for more information.
%
% See also PPFT, OptimizedPPFT, adjRadon.
%
% Yoel Shkolnisky 22/10/01

function [res1,res2]=Radon(im);

[res1,res2] = OptimizedPPFT(im);

% Inverse FFT along columns
res1 = icfftd(res1,1);
res2 = icfftd(res2,1);

% For safety - should never happen.
% No complex entries are expected in rim (the result array).
% This condition should catch programming bugs that result in complex entries in rim.
if (~isempty(find(imag(res1)>1.e-5)))
   warning(sprintf('res1 contains imaginary components of maximal value %d',max(max(imag(res1)))));
   err=1;
end
if (~isempty(find(imag(res2)>1.e-5)))
   warning(sprintf('res2 contains imaginary components of maximal value %d',max(max(imag(res2)))));
   err=1;
end

% Remove the imaginary component 00000000i from the entries of the result
% matrix. This component causes the entries of rim to be considered imaginary, 
% although the imaginary part is very close to zero.
if isreal(im)
    res1 = real(res1); 
    res2 = real(res2);
end
   
%Revision record
% 15/1/03   Yoel Shkolnisky     Use cfftd instead of column-wise cfft
