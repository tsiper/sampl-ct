function Y=FPtP(X,H)
%
% Apply the Gram operator of the pseudo-polar Fourier transform to a volume
% X. This function is called by fippft3.
%
% Yoel Shkolnisky, December 2010.

Y=applyppftfilter(X,H);