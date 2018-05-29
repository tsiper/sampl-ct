function [ theta_pp ] = ThetaPP( N , Mproj )
%THETAPP Calculates the Pseudo-Polar angles for a resolution of N
% If Mproj is not given, it will be selected as the default 2*N+2
% Mproj must be an even number

if nargin < 2
    Mproj = 2*N+2;
else
    assert(~mod(Mproj,2),'Please provide an even number for Mproj');
end

l = linspace(-N/2,N/2,Mproj/2);
theta_pp = [atan(2*l/N),fliplr(pi/2-atan(2*l/N))]*180/pi;


end

