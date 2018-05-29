% function im = adjppft3_ref(pp) 
%
% Computes the adjoint of the 3-D pseudo-polar Fourier transform.
%
% pp    4-D array of size 3x(3n+1)x(n+1)x(n+1) containing the 3-D
%       pseudo-polar Fourier transform
%
% See also ppft3_ref.
%
% Yoel Shkolnisky 26/02/03
%
% Revisions:
% Yoel Shkolnisky  19/05/2013   Renamed from adjPPFT3 to adjPPFT3_ref.


function im = adjppft3_ref(pp)

%Check if the input is of size 3x(3n+1)x(n+1)x(n+1)
n=verifyPP(pp);

%The array adjpp contains the three parts of the adjoint transform. 
%adjpp(k,:,:,:) contains the adjoint of pp(k,:,:,:)
adjpp = zeros(3,n,n,n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute the adjoint of pp(1,:,:,:)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp1 = zeros(3*n+1,n,n+1);
tmp2 = zeros(3*n+1,n,n);

% Apply adjFlipY and adjFlipZ - the adjoint of fliping along the second
% and third dimensions.
tmp  = squeeze(pp(1,:,:,:));  % Remove the singelton dimension from pp(1,:,:,:)
tmp  = flipdim(tmp,3);
tmp  = flipdim(tmp,2);


% now tmp is of size (3n+1)x(n+1)x(n+1)
% Apply adjGKN along the y direction.
for k=-3*n/2:3*n/2
    for j=-n/2:n/2
        coord = toUnaliasedCoord([k,j],[3*n+1,n]);
        U = tmp(coord{1},:,coord{2});
        tmp1(coord{1},:,coord{2}) = adjGKN3(U(:).',k);  % We use the transpose(.') since adjGKN3 requires a row vector
    end
end

% tmp1 is of size (3n+1)xnx(n+1)
% Apply adjGKN along the z direction.
for k=-3*n/2:3*n/2
    for l=-n/2:n/2-1
        coord = toUnaliasedCoord([k,l],[3*n+1,n+1]);
        V = tmp1(coord{1},coord{2},:);
        tmp2(coord{1},coord{2},:) = adjGKN3(V(:).',k); % We use the transpose(.') since adjGKN3 requires a row vector
    end
end
% now tmp2 is of size (3n+1)xnxn
% Apply adjFFT3.
tmp = prod(size(tmp2))*icfftn(tmp2);
% Apply adjEX (adjoint of padding along the X direction).
adjpp(1,:,:,:) = tmp(n+1:2*n,:,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute the adjoint of pp(2,:,:,:)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp1 = zeros(n+1,3*n+1,n+1);
tmp2 = zeros(n+1,3*n+1,n);
tmp3 = zeros(n,3*n+1,n);

% Apply adjFlipY and adjFlipZ - the adjoint of fliping along the second
% and third dimensions.
tmp  = squeeze(pp(2,:,:,:));
tmp  = flipdim(tmp,3);
tmp  = flipdim(tmp,2);

% Swap coordinates 1 and 2.
% This is the adjoint operation for the coordinate change in ppft3 for PP2.
% We perform coordinate change in ppft3 so that the pseudo-radius k will
% always be the first parameter. To compute the adjoint transform we need
% to swap the coordinates again.
% The following loop can be replaced by the line:
%     tmp1 = permute(tmp,[2 1 3]);
for k=-3*n/2:3*n/2
    for l=-n/2:n/2
        coord = toUnaliasedCoord([k,l],[3*n+1,n]);
        tmp1(coord{2},coord{1},:) = tmp(coord{1},coord{2},:);
    end
end

% now tmp1 is of size (n+1)x(3n+1)x(n+1)
% Apply adjGKN along the z direction.
for k=-3*n/2:3*n/2
    for l=-n/2:n/2
        coord = toUnaliasedCoord([k,l],[3*n+1,n]);
        V = tmp1(coord{2},coord{1},:);
        tmp2(coord{2},coord{1},:) = adjGKN3(V(:).',k);  % We use the transpose(.') since adjGKN3 requires a row vector
    end
end

% now tmp2 is of size (n+1)x(3n+1)xn
% Apply adjGKN along the x direction.
for k=-3*n/2:3*n/2
    for j=-n/2:n/2-1
        coord = toUnaliasedCoord([k,j],[3*n+1,n]);
        U = tmp2(:,coord{1},coord{2});
        tmp3(:,coord{1},coord{2}) = adjGKN3(U(:).',k); % We use the transpose(.') since adjGKN3 requires a row vector
    end
end

% now tmp3 is of size nx(3n+1)xn
% Apply adjFFT3.
tmp = prod(size(tmp3))*icfftn(tmp3);
% Apply adjEX (adjoint of padding along the X direction).
adjpp(2,:,:,:) = tmp(:,n+1:2*n,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute the adjoint of pp(3,:,:,:)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp1 = zeros(n+1,n+1,3*n+1);
tmp2 = zeros(n,n+1,3*n+1);
tmp3 = zeros(n,n,3*n+1);

% Apply adjFlipY and adjFlipZ - the adjoint of fliping along the second
% and third dimensions.
tmp  = squeeze(pp(3,:,:,:));
tmp  = flipdim(tmp,3);
tmp  = flipdim(tmp,2);

% The coordinates in pp are given as follows:
% coordinate 1 - pseudo-radius (unit steps in the z direction)
% coordinate 2 - pseudo-angle in the x direction
% coordinate 3 - pseudo-angle in the y direction.
% We rearrange the coordinates as [2 3 1] to restore the original
% coordinates oredering, where the first coordinate is x, the second is y
% and the third is z.
% The following loop can be replaced by the line:
%     tmp1 = permute(tmp,[2 3 1]);
for k=-3*n/2:3*n/2
    for l=-n/2:n/2
        coord = toUnaliasedCoord([k,l],[3*n+1,n]);
        tmp1(:,coord{2},coord{1}) = tmp(coord{1},:,coord{2});
    end
end

% now tmp1 is of size (n+1)x(n+1)(3n+1)
% Apply adjGKN along the x direction.
for k=-3*n/2:3*n/2
    for j=-n/2:n/2
        coord = toUnaliasedCoord([k,j],[3*n+1,n]);
        V = tmp1(:,coord{2},coord{1});
        tmp2(:,coord{2},coord{1}) = adjGKN3(V(:).',k); % We use the transpose(.') since adjGKN3 requires a row vector
    end
end

% now tmp2 is of size nx(n+1)x(3n+1)
% Apply adjGKN along the y direction.
for k=-3*n/2:3*n/2
    for l=-n/2:n/2-1
        coord = toUnaliasedCoord([k,l],[3*n+1,n]);
        U = tmp2(coord{2},:,coord{1});
        tmp3(coord{2},:,coord{1}) = adjGKN3(U(:).',k); % We use the transpose(.') since adjGKN3 requires a row vector
    end
end

% now tmp3 is of size nxnx(3n+1)
% Apply adjFFT3.
tmp = prod(size(tmp3))*icfftn(tmp3);
% Apply adjEX (adjoint of padding along the X direction).
adjpp(3,:,:,:) = tmp(:,:,n+1:2*n);


%%%%%%%%%%%%%%%%%%%%%%%
%Combine the 3 adjoints
%%%%%%%%%%%%%%%%%%%%%%%
im = adjpp(1,:,:,:)+adjpp(2,:,:,:)+adjpp(3,:,:,:);
im = squeeze(im);


%%%%%%%%%%%%%%%% Subroutines %%%%%%%%%%%%%%%%%

function n=verifyPP(pp)
% Verify that the input image is of size 3x(3n+1)x(n+1)x(n+1) with n even.
% If so, the function returns n. Otherwise, the function terminates with a
% proper error message.
ERR = 'pp must be of size 3x(3n+1)x(n+1)x(n+1)';
s = size(pp);

if length(s)~=4
    % The input array is not a 4-D array
    error(ERR);
end

if s(1)~=3
    % Length of first dimension is not 3
    error(ERR);
end

if (s(3)-s(4))~=0
    % The last two dimensions of pp should be n+1
    error(ERR);
end

n = s(3)-1;
if (mod(n,2)~=0)
    % n is not even
    error(ERR);
end

if s(2)~=3*n+1
    % Second dimension is not 3n+1
    error(ERR);
end