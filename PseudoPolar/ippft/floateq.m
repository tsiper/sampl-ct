function eq=floateq(x,y)
%
% Check if the real/complex numbers x and y are close up to epsilon.
% Returns 1 if equal and 0 otherwise.
%
% Yoel Shkolnisky 05/10/04
%

EPS=1.0-15;
eq= abs(x-y)<EPS;