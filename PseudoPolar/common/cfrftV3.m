% function w=cfrftV3(x,alpha)
% function w=cfrftV3(x,alpha,PQ,PZ,k)
%
% Aliased fractional Fourier transform of the sequence x.
% The FRFT is computed using O(nlogn) operations.
%
% Input:
% x       The sequence whose FRFT should be computed. Can be of odd or even
%         length. Must be a 1-D row vector.
% alpha   The parameter alpha of the fractional Fourier transform.
% PQ,PZ   Precomputed tables for avoiding unnecessary computations in
%         repeated call to the function. These tables are generated using
%         cfrftV3_precomp.
% k       Determines the value of the fractional spacing, as well as the
%         appropriate column to use in the precomputed tables. See ppft3V3.
%
% Returns the aliased FRFT with parameter alpha of the sequence x.
% The fractional Fourier transform w of the sequence x (with parameter alpha) is defined by
%                   n/2-1
%       w(j) =       sum  x(u)*exp(-2*pi*i*j*u*alpha/N),  -n/2 <= j <= n/2-1, N=length(x).
%                   u=-n/2
%
% 
% This function is the same as cfrftV2. It uses the less padding (3m as in
% the paper) and therefore it is more memory efficient. This may cause the
% lengths of the sequences to be non-optimal for Matlab's FFT function. The
% function cfrftV2 uses more memory (for padding) but uses FFTs of dyadic
% length. 
%
% If five parameters are provided, then alpha is ignored. Just pass,
% for example, zero.
%
% Yoel Shkolnisky 22/10/01
%
% Revisions:
% Yoel Shkolnisky 21/05/2013 Chage the code to work with column vectors.
%          Allow x to be a matrix, in which case the transform is applied
%          on all columns.

function w=cfrftV3(x,alpha,PQ,PZ,k)

m=length(x);
lm= -fix(m/2);
hm= fix((m-0.5)/2);
ofs=floor(3*m/2)+1;

if nargin==5
    % load weights from the precomputed structure
    q=(PQ(:,k));
    Z=(PZ(:,k));
else
    % computed required weights
    j=lm:hm;
    j2= -m:m;
    E=1i*pi*alpha;
    q=exp(-E*j.^2/m);
    q=q(:);
    
    z=zeros(3*m,1);
    l=-m+ofs;
    z(l:l+length(j2)-1)=exp(E*j2.^2/m);
    Z=fft(z);
end

sz2=size(x,2);
y=bsxfun(@times,x,q);
y=[zeros(m,sz2);y;zeros(m,sz2)];
Y=fft(y);
W=bsxfun(@times,Y,Z);
w=ifft(W);
w=w([ofs:3*m 1:ofs-1],:);
w=w(lm+ofs:hm+ofs,:);
w=bsxfun(@times,w,q);
