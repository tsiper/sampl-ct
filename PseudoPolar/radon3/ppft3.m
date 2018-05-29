% function pp=ppft3(im);
%
% Fast algorithm for computing the 3-D pseudo-polar Fourier transform.
% The computation requires O(n^3logn) operations.
%
% The function computes the 3-D pseudo-polar Fourier transform according to
% the algorithm given in
% "A. Averbuch and Y. Shkolnisky. 3D Fourier based discrete Radon
% transform. Applied and Computational Harmonic Analysis, 15(1):33-69,
% 2003."
%
% This implementation follows exactly the pseudo-code and notations given
% in the paper. 
%
% Input:
%     im    3-D image of size nxnxn (n even).
%           First  index - x direction
%           Second index - y direction
%           Third  index - z direction
%     The indices in each direction are assumed to be from -n/2 to n/2-1 and
%     not from 1 to n.
%
% Output:
%     pp - 4-D array of size 3x(3n+1)x(n+1)x(n+1) containing the 3-D
%     pseudo-polar Fourier transform.
%     The array pp contains the following Fourier samples:
%        pp(1,k,l,j) = FI(k,-2lk/n,-2jk/n)
%        pp(2,k,l,j) = FI(-2lk/n,k,-2jk/n)
%        pp(3,k,l,j) = FI(-2lk/n,-2jk/n,k)
%        where
%             l,j = -n/2,...,n/2   k=-3n/2,...3n/2
%        and
%                           n/2-1   n/2-1   n/2-1
%             FI(ox,oy,oz)=  sum     sum     sum  I(u,v,w)exp(-2*pi*i(u*ox+v*oy+w*oz)/m)   m=3n+1
%                           u=-n/2  v=-n/2  w=-n/2
%
% See also ppft3_ref.
%
% Yoel Shkolnisky 30/01/03
%
% Revisions:
% Yoel Shkolnisky   19/05/2013    Renamed from OptimizedPPFT3 to PPFT3.
% Yoel Shkolnisky   21/05/2013    Vectorize for speed (about factor 4).

function pp=ppft3(im)

% verify that the input is a 3D image of size nxnxn
verifyImage(im);

% Initialize output data structure
s=size(im);
n=s(1); % at this point n is even
m = 3*n+1;
pp  = zeros(3,3*n+1,n+1,n+1);
tmp = zeros(3*n+1,n+1,n+1);
alpha = 2*(n+1)/(n*m);

[PQ,PZ]=cfrftV3_precomp(n+1,alpha.*(-3*n/2:3*n/2));

% Compute the pseudo-polar Fourier transform PP1

pim  = cat(1,zeros(n,n,n),im,zeros(n+1,n,n)); % pad the image to size m along the x direction
fim  = cfftd(pim,1);
tmp1 = zeros(m,n,n+1); % intermediate result after the first resampling. Referred as T1 in the paper.

ofs_m=floor(m/2)+1;
%ofs_n=floor(n/2)+1;

for k=-3*n/2:3*n/2
%     for l=-n/2:n/2-1
%         U = fim(k+ofs_m,l+ofs_n,:);
%         tmp1(k+ofs_m,l+ofs_n,:) = cfrftV3([U(:); 0],0,PQ,PZ,k+3*n/2+1);
%     end
    U=squeeze(fim(k+ofs_m,:,:));
    U=U.';
    F=cfrftV3([U; zeros(1,n)],0,PQ,PZ,k+3*n/2+1);
    tmp1(k+ofs_m,:,:)=F.';
end

for k=-3*n/2:3*n/2
%     for j=-n/2:n/2
%         V = tmp1(k+ofs_m,:,j+ofs_n);
%         tmp(k+ofs_m,:,j+ofs_n) = cfrftV3([V(:); 0],0,PQ,PZ,k+3*n/2+1);
%     end    
    V=squeeze(tmp1(k+ofs_m,:,:)); 
    F=cfrftV3([V; zeros(1,n+1)],0,PQ,PZ,k+3*n/2+1);
    tmp(k+ofs_m,:,:)=F;
end

tmp = flipdim(tmp,2);
pp(1,:,:,:) = flipdim(tmp,3);

% Compute the pseudo-polar Fourier transform PP2

pim  = cat(2,zeros(n,n,n),im,zeros(n,n+1,n)); % pad the image to size m along the y direction
fim  = cfftd(pim,2);
tmp1 = zeros(n+1,m,n); % intermediate result after the first resampling. Referred as T2 in the paper.

% The loop order (k,l,j) differs from PP1 and PP3 to keep consistency with
% the paper.
for k=-3*n/2:3*n/2
%     for j=-n/2:n/2-1
%         U = fim(:,k+ofs_m,j+ofs_n);
%         tmp1(:,k+ofs_m,j+ofs_n) = cfrftV3([U(:); 0],0,PQ,PZ,k+3*n/2+1);
%     end
    U=squeeze(fim(:,k+ofs_m,:));
    F=cfrftV3([U; zeros(1,n)],0,PQ,PZ,k+3*n/2+1);
    tmp1(:,k+ofs_m,:)=F;
end

for k=-3*n/2:3*n/2
%     for l=-n/2:n/2
%         V = tmp1(l+ofs_n,k+ofs_m,:);
%         tmp(k+ofs_m,l+ofs_n,:) = cfrftV3([V(:); 0],0,PQ,PZ,k+3*n/2+1);
%     end
    V=squeeze(tmp1(:,k+ofs_m,:));
    V=V.';
    F=cfrftV3([V; zeros(1,n+1)],0,PQ,PZ,k+3*n/2+1);
    tmp(k+ofs_m,:,:)=F.';

end

tmp = flipdim(tmp,2);
pp(2,:,:,:) = flipdim(tmp,3);

% Compute the pseudo-polar Fourier transform PP3

pim  = cat(3,zeros(n,n,n),im,zeros(n,n,n+1)); % pad the image to size m along the z direction
fim  = cfftd(pim,3);
tmp1 = zeros(n,n+1,m); % intermediate result after the first resampling. Referred as T3 in the paper.

for k=-3*n/2:3*n/2
%     for l=-n/2:n/2-1
%         U = fim(l+ofs_n,:,k+ofs_m);
%         tmp1(l+ofs_n,:,k+ofs_m) = cfrftV3([U(:); 0],0,PQ,PZ,k+3*n/2+1);
%     end
    U=squeeze(fim(:,:,k+ofs_m));
    U=U.';
    F=cfrftV3([U; zeros(1,n)],0,PQ,PZ,k+3*n/2+1);
    tmp1(:,:,k+ofs_m)=F.';    
end

for k=-3*n/2:3*n/2
%     for j=-n/2:n/2
%         V = tmp1(:,j+ofs_n,k+ofs_m);
%         tmp(k+ofs_m,:,j+ofs_n) = cfrftV3([V(:); 0],0,PQ,PZ,k+3*n/2+1);
%     end
    V=squeeze(tmp1(:,:,k+ofs_m));    
    F=cfrftV3([V; zeros(1,n+1)],0,PQ,PZ,k+3*n/2+1);
    tmp(k+ofs_m,:,:)=F;    
end

tmp = flipdim(tmp,2);
pp(3,:,:,:) = flipdim(tmp,3);