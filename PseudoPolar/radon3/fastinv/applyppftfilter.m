function Y=applyppftfilter(X,H)
%
% Apply the convolution filter that corresponds to the gram operator of the
% pseudo-polar Fourier transform to the volume X.
%
% X should be a volume of size nxnxn, and H should be a convolution filter
% for size n, as generared by makeppftfilter(n).
%
% Yoel Shkolnisky, December 2010.

verifyImage(X);
n=size(X,1);
Xp=zeros(3*n,3*n,3*n);
Xp(n+1:2*n,n+1:2*n,n+1:2*n)=X;

hp=zeros(3*n,3*n,3*n);
hp(n/2+1:5*n/2,n/2+1:5*n/2,n/2+1:5*n/2)=H.filter;
Hhat=fftn(ifftshift(hp));

Y=fftshift(ifftn(fftn(ifftshift(Xp)).*Hhat));
Y=Y(n+1:2*n,n+1:2*n,n+1:2*n);
