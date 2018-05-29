function H=makeppftfilter(n,precond,precision)
%
% Compute the convolution kernel required for applying the Gram operator of
% the pseudo-polar Fourier transform. The filter is created in the given
% precision ('single' or 'double').
%
% nxnxn are the dimensions on which the ppft is applied.
%
% Set precond to nonzero to compute the convolution filter of the
% preconditioned gram operator.
%
% The computed filter H is an object with the following fields:
%   filter     The data of the filter
%   precond    1 if the filter is for the preconditioned gram operator and
%              0 otherwise.
%   precision  'single' or 'double'. Default is 'double'.
%
% Yoel Shkolnisky, December 2010.
%
% Revisions:
% April 2014  Y.S.  Add precision variable.


if ~exist('precond','var')
    precond=0;
end

if ~exist('precision','var')
    precision='double';
end

m=3*n+1;

alpha=ones(3*n+1,n+1,n+1);
c=2*(n+1)/(m*n);
if precond
    for k=-3*n/2:3*n/2
        if k==0
            alpha(k+3*n/2+1,:,:)=1/m^4;
        else
            alpha(k+3*n/2+1,:,:)=(abs(k*c)).^2;
        end
    end
end

% Generate the frequencies of the pseudo-polar Fourier transform that
% correspond to the first sector of the PP grid (PP(1,:,:,)). 
% The other two sectors are obtained by shuffling these frequencies, as
% done
[ox,oy,oz]=ndgrid(-3*n/2:3*n/2, -n/2:n/2, -n/2:n/2); 
omega=[ox(:) oy(:) oz(:)];
omega(:,2)=omega(:,2).*(-2).*omega(:,1)./n;
omega(:,3)=omega(:,3).*(-2).*omega(:,1)./n;

% Compute h1

h1b=nufft_3d(alpha,2.*n.*omega/m,precision,2*n);

% Compute h2
omega=[omega(:,2) omega(:,1) omega(:,3)];
h2b=nufft_3d(alpha,2.*n.*omega/m,precision,2*n);

% Compute h3
omega=[omega(:,1) omega(:,3) omega(:,2)];
h3b=nufft_3d(alpha,2.*n.*omega/m,precision,2*n);

h=h1b+h2b+h3b;

if strcmpi(precision,'single')
    h=single(h);
end

H.filter=h;
H.precond=precond;


% The following version works for any n, but is not needed since PPFT
% requires n to be even.
% % precision='double';
% % m=3*n+1;
% % alpha=ones(3*n+1,n+1,n+1);
% % 
% % % Generate the frequencies of the pseudo-polar Fourier transform that
% % % correspond to the first sector of the PP grid (PP(1,:,:,)). 
% % % The other two sectors are obtained by shuffling these frequencies, as
% % % done below.
% % [ox,oy,oz]=ndgrid(-fix(3*n/2):fix(3*n/2),-fix(n/2):fix(n/2),-fix(n/2):fix(n/2)); 
% % omega=[ox(:) oy(:) oz(:)];
% % omega(:,2)=omega(:,2).*(-2).*omega(:,1)./n;
% % omega(:,3)=omega(:,3).*(-2).*omega(:,1)./n;
% % 
% % % Compute h1
% % h1b=nufft_3d(alpha,2.*n.*omega/m,precision,2*n);
% % 
% % % Compute h2
% % omega=[omega(:,2) omega(:,1) omega(:,3)];
% % h2b=nufft_3d(alpha,2.*n.*omega/m,precision,2*n);
% % 
% % % Compute h3
% % omega=[omega(:,1) omega(:,3) omega(:,2)];
% % h3b=nufft_3d(alpha,2.*n.*omega/m,precision,2*n);
% % 
% % % Compute the Fourier transform of the convolution filter.
% % h=h1b+h2b+h3b;
% % hp=zeros(3*n,3*n,3*n);
% % %hp(n/2+1:5*n/2,n/2+1:5*n/2,n/2+1:5*n/2)=h;
% % idxl=toUnaliasedIdx(-n,3*n);
% % idxh=toUnaliasedIdx(n,3*n)-1;
% % hp(idxl:idxh,idxl:idxh,idxl:idxh)=h;
% % H=fftn(ifftshift(hp));


