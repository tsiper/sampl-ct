% function [sec1,sec2] = extractPPsectors(ppim)
%
% Extract a combined pseudo-polar image into 2 pseudo-polar sectors.
% ppim must be of size (2n+1)x(2n+2).
%
% See combinePPsectors for description of sec1, sec2 and ppim.
%
% Yoel Shkolnisky 22/10/01


function [sec1,sec2] = extractPPsectors(ppim)

%Check if the input is of valid size (2n+1)x(2n+2)
s = size(ppim);
if (mod(s(1),2)~=1) | (mod(s(2),2)~=0)
   error('Input parameter must be of size (2n+1)x(2n+2)');
end
n = [(s(1)-1)/2;(s(2)-2)/2];
if (n(1)~=n(2))
   error('Input parameter must be of size (2n+1)x(2n+2)');
end
n = n(1);

sec1 = ppim(:,1:n+1);
sec2 = fliplr(ppim(:,n+2:2*n+2));


