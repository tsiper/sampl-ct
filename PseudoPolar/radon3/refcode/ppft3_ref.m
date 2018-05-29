% function pp = ppft3_ref(im)
%
% Fast algorithm for computing the 3-D pseudo-polar Fourier transform.
% The computation requires O(n^3logn) operations.
%
% The function computes the 3-D pseudo-polar Fourier transform according to
% the algorithm given in
% "A. Averbuch and Y. Shkolnisky. 3D Fourier based discrete Radon
% transform. Applied and Computational Harmonic Analysis, 15(1):33-69,
% 2003."
% This implementation follows exactly the pseudo-code and notations given
% in the paper. See ppft3 for an optimized version of the algorithm.
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
%
% Yoel Shkolnisky 30/01/03
%
% Revisions:
% Yoel Shkolnisky 19/5/2013  Renamed from ppft3 to ppft3_ref.

function pp=ppft3_ref(im)

% verify that the input is a 3D image of size nxnxn
verifyImage(im);

% Initialize output data structure
s=size(im);
n=s(1); % at this point n is even
m = 3*n+1;
pp  = zeros(3,3*n+1,n+1,n+1);
tmp = zeros(3*n+1,n+1,n+1);

% Compute the pseudo-polar Fourier transform PP1

pim  = cat(1,zeros(n,n,n),im,zeros(n+1,n,n)); % pad the image to size m along the x direction
fim  = cfftn(pim);
tmp1 = zeros(m,n,n+1); % intermediate result after the first resampling. Referred as T1 in the paper.

for k=-3*n/2:3*n/2
    for l=-n/2:n/2-1
        coord = toUnaliasedCoord([k,l],[m,n]);
        U = fim(coord{1},coord{2},:);
        tmp1(coord{1},coord{2},:) = GKN3(U(:).',k); % We use the transpose(.') since GKN3 requires a row vector
    end
end

for k=-3*n/2:3*n/2
    for j=-n/2:n/2
        coord = toUnaliasedCoord([k,j],[m,n]);
        V = tmp1(coord{1},:,coord{2});
        tmp(coord{1},:,coord{2}) = GKN3(V(:).',k);  % We use the transpose(.') since GKN3 requires a row vector
    end
end

tmp = flipdim(tmp,2);
pp(1,:,:,:) = flipdim(tmp,3);

% Compute the pseudo-polar Fourier transform PP2

pim  = cat(2,zeros(n,n,n),im,zeros(n,n+1,n)); % pad the image to size m along the y direction
fim  = cfftn(pim);
tmp1 = zeros(n+1,m,n); % intermediate result after the first resampling. Referred as T2 in the paper.

% The loop order (k,l,j) differs from PP1 and PP3 to keep consistency with
% the paper.
for k=-3*n/2:3*n/2
    for j=-n/2:n/2-1
        coord = toUnaliasedCoord([k,j],[m,n]);
        U = fim(:,coord{1},coord{2});
        tmp1(:,coord{1},coord{2}) = GKN3(U(:).',k); % We use the transpose(.') since GKN3 requires a row vector
    end
end

for k=-3*n/2:3*n/2
    for l=-n/2:n/2
        coord = toUnaliasedCoord([k,l],[m,n]);
        V = tmp1(coord{2},coord{1},:);
        tmp(coord{1},coord{2},:) = GKN3(V(:).',k);  % We use the transpose(.') since GKN3 requires a row vector
    end
end

tmp = flipdim(tmp,2);
pp(2,:,:,:) = flipdim(tmp,3);

% Compute the pseudo-polar Fourier transform PP3

pim  = cat(3,zeros(n,n,n),im,zeros(n,n,n+1)); % pad the image to size m along the z direction
fim  = cfftn(pim);
tmp1 = zeros(n,n+1,m); % intermediate result after the first resampling. Referred as T3 in the paper.

for k=-3*n/2:3*n/2
    for l=-n/2:n/2-1
        coord = toUnaliasedCoord([k,l],[m,n]);
        U = fim(coord{2},:,coord{1});
        tmp1(coord{2},:,coord{1}) = GKN3(U(:).',k); % We use the transpose(.') since GKN3 requires a row vector
    end
end

for k=-3*n/2:3*n/2
    for j=-n/2:n/2
        coord = toUnaliasedCoord([k,j],[m,n]);
        V = tmp1(:,coord{2},coord{1});
        tmp(coord{1},:,coord{2}) = GKN3(V(:).',k);  % We use the transpose(.') since GKN3 requires a row vector
    end
end

tmp = flipdim(tmp,2);
pp(3,:,:,:) = flipdim(tmp,3);
