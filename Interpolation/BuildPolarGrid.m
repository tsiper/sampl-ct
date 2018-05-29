function [ PolarGrid , theta ] = BuildPolarGrid( Nray,Mproj_Theta )
%BUILDPPMTX This function builds the polar grid, for computing the from parallel
%beam to a cartesian grid. We have Nray - detectors and Mproj projection angles
%spanning from 0 to 180 degrees
% The input theta is in degrees

% PPgrid = zeros(m,n+1);
% PPgrid = PPgrid;

r     = -Nray/2+1/2:Nray/2-1/2;

% If Mproj_Theta is a scalar we calculate the angles, otherwise we take it as a
% vector of theta
if isscalar(Mproj_Theta)
    Mproj = Mproj_Theta;
    theta = (0:1/Mproj:1-1/Mproj)*pi;
else
    theta = deg2rad(Mproj_Theta);
end

[THETA,R] = meshgrid(theta,r);

% R(:,end/2+1:end) = R(:,end/2+1:end)+1;

R = fliplr(R);

X = -R.*sin(THETA);
Y = R.*cos(THETA);

PolarGrid = [Y(:),X(:)];

end

