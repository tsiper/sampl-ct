function ippftprecomp(n)
%
% Compute the tables used in the inversion of the pseudo-polar of an image
% of size nxn.
% The function generates two mat files named c_[n] and d_[n] (for example
% c_256 and d_256.
%
% Yoel Shkolnisky 11/10/04

[c,d]=ippftconsts(n);
fname=sprintf('./d_%d',n);
save(fname,'d');
fname=sprintf('./c_%d',n);
save(fname,'c');

