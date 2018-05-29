function Y=applyppfilter_gpu(gX,gH,precision)
%
% GPU version of applyppfilter.m.
%
% Yoel Shkolnisky, 21/05/2013

n=size(gX,1);
gXp=gpuArray.zeros(3*n,3*n,3*n,precision);
gXp(n+1:2*n,n+1:2*n,n+1:2*n)=gX;
Y=fftshift(ifftn(fftn(ifftshift(gXp)).*gH));
Y=Y(n+1:2*n,n+1:2*n,n+1:2*n);

