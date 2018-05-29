function [ R ] = ppRadon( I, N)
%PP_RADON We start with the Image I and use equi-slope projections to produce
%the singoram with N angles the image I is supposed to be a perfectly square
%image

[n,m] = size(I);

if (n ~= m) || (mod(n,2) ~= 0)
    error('Please insert a square image I~(NxN), where N is even');
end
M = 2*n;
% I = @(u,v) Image(u+n/2+1,v+n/2+1);

% The linear slope parameters y=sx+t / x=sy+t
% t = -n:n-1;
% s = linspace(-1,1,N);
% s = 0;

% R1 = zeros(length(s),length(t));
% R2 = R1;
% 
% Ix = padarray(I,n/2);
% Iy = padarray(I,[0,n/2]);
% 
% figure;
% subplot(211);imagesc(Ix);
% subplot(212);imagesc(Iy);

% computing Radon transform for y and x
% for i =1:length(s)
%     for j=1:length(t)
%         for u=-n/2:n/2-1
%             for v=-n/2:n/2-1
%                 R1(i,j) = R1(i,j) + Ix(u+n+1,v+n/2+1)*diric(2*pi/M*(s(i)*u+t(j)-v),M);
%                 R2(i,j) = R2(i,j) + Iy(u+n/2+1,v+n+1)*diric(2*pi/M*(s(i)*v+t(j)-u),M);
% %                 R3(i,j) = R3(i,j) + I(u+n/2+1,v+n/2+1)*diric(2*pi*(s(i)*u+t(j)-v),m);
% %                 R4(i,j) = R4(i,j) + I(u+n/2+1,v+n/2+1)*diric(2*pi*(s(i)*v+t(j)-u),m);
%             end
%         end
%     end
% end
% 
% R = [R1,R2];

% % Equispaced angled
% thetas = (0:1/N:1-1/N)*180;

% Equisloped angles
% thetas = [atan(-1:2/N:0-2/N),atan(0:2/N:1-2/N),mod(acot(1:-2/N:-1+2/N),pi)+pi];
thetas = [atan(-1:2/N:1-2/N),atan(-1:2/N:1-2/N)+pi/2];
thetas = thetas/pi*180;

R = zeros(M,N);
for j=1:2*N
    R(:,j) = Rt(I,thetas(j));
end


end


