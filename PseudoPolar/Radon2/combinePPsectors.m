% function ppim = combinePPsectors(sec1,sec2)
%
% Take the two sectors of samples on the pseudo-polar grid and combine them into
% a single pseudo-polar image.
% 
% sec1, sec2    Two matrices of samples on the pseudo-polar grid. sec1(k,l),sec2(k,l) 
%               represent the samples at radius "k" and angle "l".
%
% The angles along sec1 columns correspond to theta=arctan(2l/n) l=-n/2...n/2 (from -pi/4 to pi/4).
% The angles along sec2 columns correspond to theta=pi/2-arctan(2l/n) l=-n/2...n/2 (from 3p/4 to pi/4).
%
% For the combined image ppim(k,l), k is radius and l is angle.
% The angles in the combined matrix ppim are from -pi/4 to 3pi/4. More specifically,
% the angles are ordered  as -pi/4...pi/4,pi/4...3pi/4.
% 
% We do not adjust the meaning of the parameter k (pseudo-radius) for each of the sectors. 
% For -pi/4...pi/4 it measures the intercept along the y-axis. For pi/4...3pi/4 it measures 
% the intercept along the x-axis.
%
% Yoel Shkolnisky 22/10/01

function ppim = combinePPsectors(sec1,sec2);

s1 = size(sec1);
s2 = size(sec2);

if (length(s1)~=2 | length(s2)~=2)
   error('sec1 and sec2 must be 2D matrices');
end
if (s1(1)~=s2(1) | s1(2)~=s2(2))
   error('sec1 and sec2 must have the same size');
end

a = sec1;
b = fliplr(sec2);

ppim = [a b];
 
