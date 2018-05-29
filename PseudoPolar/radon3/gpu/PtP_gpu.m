function Y=PtP_gpu(gX,gH,precision)
%
% Apply the Gram operator of the pseudo-polar Fourier transform to a volume
% X. This function is called by fippft3_gpu.
%
% Yoel Shkolnisky, 21/05/2013

Y=applyppfilter_gpu(gX,gH,precision);