function numOfScales = defaultNumberOfScales(l)
%DEFAULTNUMBEROFSCALES compute default number of scales
%
% INPUT:
%  l	            (vector) size of image
%
% OUTPUT:
%  numOfScales      (int) number of scales
%
%--------------------------------------------------------------------------
% Sören Häuser ~ FFST ~ 2014-07-22 ~ last edited: 2014-07-22 (Sören Häuser)

    numOfScales = floor(0.5 * log2(max(l)));
end