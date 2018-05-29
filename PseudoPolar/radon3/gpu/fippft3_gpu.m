function [Y,flag,residual,iter] = fippft3_gpu(pp,ErrTol,MaxIts,verbose,precision)
%
% Fast inverse pseudo-polar Fourier transform using GPU acceleration.
% See fippft3.m for more information
%
% Yoel Shkolnisky, 21/05/2013.
%

if nargin<5
   precision='double';
end

if nargin<4
   verbose = 0;
end

if nargin<3
   MaxIts = 10;
end

if nargin<2
   ErrTol = 1.e-2;
end

n=verifyPP(pp);
filtname=getppfiltname(n,precision);
H=loadppftfilter(filtname);
if ~isstruct(H)
    fprintf('Filter not found. Generating filter\n');
    H=makeppftfilter(n,1,precision);
    saveppftfilter(H,filtname);
end

if verbose
    fprintf('GPU PRECISION=%s\n',precision);
    if strcmpi(precision,'single')
        fprintf('Using single precision GPU calculations!!!\n');
    end
end


% Take the FFT of the filter for fast convolution in PtP_gpu.
hp=zeros(3*n,3*n,3*n);
hp(n/2+1:5*n/2,n/2+1:5*n/2,n/2+1:5*n/2)=H.filter;
Hhat=fftn(ifftshift(hp));

if strcmpi(precision,'single')
    gHhat=gpuArray(single(Hhat));
else
    gHhat=gpuArray(double(Hhat));
end

temp = precondadjppft3(pp);

if strcmpi(precision,'single')
    gtemp=gpuArray(single(temp));
else
    gtemp=gpuArray(double(temp));
end
[gY,flag,residual,iter] = CG('PtP_gpu',gtemp,{gHhat,precision},ErrTol,MaxIts,zeros(size(temp)),verbose);
Y=gather(gY);
if flag
   warning ('Inversion did not converge. Residual error %-2.5e',residual);
end
