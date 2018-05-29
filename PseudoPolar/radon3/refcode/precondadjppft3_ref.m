% function im = precondadjppft3_ref(pp) 
%
% Computes the preconditioned adjoint of the 3-D pseudo-polar Fourier transform.
%
% pp    4-D array of size 3x(3n+1)x(n+1)x(n+1) containing the 3-D
%       pseudo-polar Fourier transform
%
% See also ppft3, adjppft3_ref, adjppft3.
%
% Yoel Shkolnisky 28/02/03
%
% Revision:
% Yoel Shkolnisky  02/07/2007 
%    n was replaced with n+1. This does not affect anything since n is even
%    (since we take 3*n/2 and expect to get an integer), and so
%    floor((n+1)/2) is exactly the same as floor(n/2).
% Yoel Shkolnisky  19/05/2013
%     Rename precondAdjPPFT3 to precondAdjPPFT3_ref.


function im = precondadjppft3_ref(pp)

%Check if the input is of size 3x(3n+1)x(n+1)x(n+1)
n=verifyPP(pp);

%The array adjpp contains the three parts of the adjoint transform. 
%adjpp(k,:,:,:) contains the adjoint of pp(k,:,:,:)
adjpp = zeros(3,n,n,n);
m=3*n+1;
alpha = 2*(n+1)/(n*m);

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
% Apply adjoint of FRFT along the y direction.
for k=-3*n/2:3*n/2
    for j=-n/2:n/2
        coord = toUnaliasedCoord([k,j],[3*n+1,n+1]);
        U = tmp(coord{1},:,coord{2});
        F = mult(k,alpha,n).*cfrft(U(:).',-k*alpha);  % (:).' is used since cfrft requries a row vector
        tmp1(coord{1},:,coord{2}) = F(1:n); % Truncate the last element of F
    end
end

% tmp1 is of size (3n+1)xnx(n+1)
% Apply adjoint of FRFT along the z direction.
for k=-3*n/2:3*n/2
    for l=-n/2:n/2-1
        coord = toUnaliasedCoord([k,l],[3*n+1,n+1]);
        V = tmp1(coord{1},coord{2},:);
        F = mult(k,alpha,n).*cfrft(V(:).',-k*alpha);
        tmp2(coord{1},coord{2},:) = F(1:n);
    end
end
% now tmp2 is of size (3n+1)xnxn
% Apply adjFFT along the x direction.
tmp = m*icfftd(tmp2,1);
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
tmp1 = permute(tmp,[2 1 3]);

% now tmp1 is of size (n+1)x(3n+1)x(n+1)
% Apply the adjoint of FRFT along the z direction.
for k=-3*n/2:3*n/2
    for l=-n/2:n/2
        coord = toUnaliasedCoord([k,l],[3*n+1,n+1]);
        V = tmp1(coord{2},coord{1},:);
        F = mult(k,alpha,n).*cfrft(V(:).',-k*alpha);
        tmp2(coord{2},coord{1},:) = F(1:n);
    end
end

% now tmp2 is of size (n+1)x(3n+1)xn
% Apply the adjoint of FRFT along the x direction.
for k=-3*n/2:3*n/2
    for j=-n/2:n/2-1
        coord = toUnaliasedCoord([k,j],[3*n+1,n+1]);
        U = tmp2(:,coord{1},coord{2});
        F = mult(k,alpha,n).*cfrft(U(:).',-k*alpha);
        tmp3(:,coord{1},coord{2}) = F(1:n);
    end
end

% now tmp3 is of size nx(3n+1)xn
% Apply adjFFT along the y direction.
tmp = m*icfftd(tmp3,2);
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
tmp1 = permute(tmp,[2 3 1]);

% now tmp1 is of size (n+1)x(n+1)(3n+1)
% Apply the adjoint of FRFT along the x direction.
for k=-3*n/2:3*n/2
    for j=-n/2:n/2
        coord = toUnaliasedCoord([k,j],[3*n+1,n+1]);
        V = tmp1(:,coord{2},coord{1});
        F = mult(k,alpha,n).*cfrft(V(:).',-k*alpha);
        tmp2(:,coord{2},coord{1}) = F(1:n);
    end
end

% now tmp2 is of size nx(n+1)x(3n+1)
% Apply the adjoint of FRFT along the y direction.
for k=-3*n/2:3*n/2
    for l=-n/2:n/2-1
        coord = toUnaliasedCoord([k,l],[3*n+1,n+1]);
        U = tmp2(coord{2},:,coord{1});
        F = mult(k,alpha,n).*cfrft(U(:).',-k*alpha);
        tmp3(coord{2},:,coord{1}) = F(1:n);
    end
end

% now tmp3 is of size nxnx(3n+1)
% Apply adjFFT along.
tmp = m*icfftd(tmp3,3);
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
% If so, the function return n. Otherwise, the function terminates with
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


function v = mult(k,alpha,n)
% Compute preconditioning factor for row k
if k==0
    v = 1/((3*n+1)^2);
else
    v = abs(k*alpha);
end
   
