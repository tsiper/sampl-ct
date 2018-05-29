function [ img, simfid ] = squares_fast( N,sqrmat,kvec,percent,noplot )

%  Produce sythetic data of squares phantom.
%  Simulate of both image and fids.
%  Compute both simulated fid and fourier transformed fid.
%  For computation of synthetic FID see routine 'sqr2'.
%
%  Input:
%
%    N		- number of points (power of 2)
%    sqrmat(6,M)	- M squares. Each square has 6 parameters:
% 			sqrmat(0,*) - amplitude
% 			sqrmat(1,*) - x-length (integer)
% 			sqrmat(2,*) - y-length (integer)
% 			sqrmat(3,*) - x-shift  (integer)
% 			sqrmat(4,*) - y-shift  (integer)
% 			sqrmat(5,*) - rotation angle  (degrees)
%
%    kvec(2,N1)	- 'simfid' is sampled at
% 		  these points of the k-space. For non-liniar sampling.
%
%  Keywords:
%
%    \percent	- (optional) when exists: lengths and shifts are given
% 		  as a percent of total number of points (N).
%    noplot	- 1 = don't plot (optional, default = 0)
%
%  Output:
%
%    img(N,N)	- image domain. squares.
%    simfid(N1)	- simulated fid (k-space). Computed by including an appropriate
% 		  2D-sinc function for each square.

temp_sqrmat = sqrmat;
sqrmat(2,:) = temp_sqrmat(3,:);
sqrmat(3,:) = temp_sqrmat(2,:);
sqrmat(4,:) = -temp_sqrmat(5,:);
sqrmat(5,:) = -temp_sqrmat(4,:);


sz = size(sqrmat);
M = sz(2);
N1 = length(kvec);
N2 = N/2;

img = zeros(N);
simfid = zeros(1,N1);
[X,Y] = meshgrid(2:1:N+1, 2:1:N+1);

h = waitbar(0,'Generating Phantom...');
for sqri = 1:M
    tmp_image = zeros(N);
    amp = sqrmat(1,sqri);
    if percent == 1
        la = sqrmat(2,sqri).*N;
        lb = sqrmat(3,sqri).*N;
        shiftx = sqrmat(4,sqri).*N;
        shifty = sqrmat(5,sqri).*N;
    else
        la = sqrmat(2,sqri);
        lb = sqrmat(3,sqri);
        shiftx = sqrmat(4,sqri);
        shifty = sqrmat(5,sqri);
    end
    ang = sqrmat(6,sqri);
    angrad = ang*pi/180;
    rotmat = [cos(angrad) -sin(angrad) ; sin(angrad), cos(angrad)]';
    rotmat1 = rotmat';
    
    rlena=2*la+1;
    rlenb=2*lb+1;
    V =rotmat1*[X(:).'-N2-shiftx-1 ; Y(:).'-N2-shifty-1];
    Vx = reshape(V(1,:),N,N);
    Vy = reshape(V(2,:),N,N);
    tmp_image = amp.*((abs(Vx) <= la) & (abs(Vy) <= lb));
    img = img + tmp_image.';
    krot = rotmat1*kvec;
    tmpfid = amp.*sinc(krot(1,:).*rlena./N)./sinc(krot(1,:)./N).*sinc(krot(2,:).*rlenb./N)./sinc(krot(2,:)./N).*rlena.*rlenb.*exp(-2*1i*pi*kvec(1,:)*shiftx./N).*exp(-2*1i*pi*kvec(2,:)*shifty/N);
    %     tmpfid = amp.*sinc(krot(1,:).*rlena./N).*sinc(krot(2,:).*rlenb./N).*rlena.*rlenb.*exp(-2*1i*pi*kvec(1,:)*shiftx./N).*exp(-2*1i*pi*kvec(2,:)*shifty/N);
    simfid = simfid + tmpfid;
    waitbar(sqri / M)
end
close(h)
end

