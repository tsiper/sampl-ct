function [ PPGrid , theta_vec ] = BuildPPGrid( n )
%BUILDPPMTX This function builds the PP grid, for computing the pseudo-polar
%transform, according to equation 2 from Shkolnitzky (I) paper:

% The dimensions of the measurements
% close all;
% n =5;

m = 2*n+1;

% PPgrid = zeros(m,n+1);
% PPgrid = PPgrid;

l = -n/2:n/2;
k = -n:n;

[X,Y] = meshgrid(l,k);


Z = -2*X.*Y/n;


PPGrid = [Y(:),Z(:)];
PPGrid = [PPGrid; Z(:),Y(:)];

% 
% % Deubg Plot
% for j = 0:(2*n+2)-1
%     scatter(PPgrid((j*m+1):(j+1)*m,1),PPgrid((j*m+1):(j+1)*m,2));
%     xlim([-n-1,n+1]);
%     ylim([-n-1,n+1]);
%     pause(0.1);
%     hold on;
% end

% PPgrid = [PPgrid, atand(PPgrid(:,2)./(PPgrid(:,1)+eps))];

% figure


%% Fixing the angles so that comply with the PP regime
PPGridFix = PPGrid;
for j=1:n+1
    range     = (j-1)*m+1 : j*m;
    new_range = (n-j+1)*m+1 : (n-j+2)*m;
    PPGridFix( range , : ) = PPGrid( new_range , : );
end
% end
% 

% Now flipping the last quadrant
% for j=3*n/2+2 : 2*n+2
%     range     = (j-1)*m+1 : j*m;
%     PPGridFix(range,:)
PPGridFix((3*n/2+1)*m+1:end,:) = (PPGrid((3*n/2+1)*m+1:end,:));
% end

PPGrid = PPGridFix;
j = 1:m:m*(2*n+2);
theta_vec = atand(PPGrid(j,2)./(PPGrid(j,1)+eps));

% 