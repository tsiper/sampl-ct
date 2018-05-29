function [ GridCoordinates ] = ConstructGridRadial0( Nbins, NSpokes)

[rr, ss] = meshgrid(0:Nbins-1,0:NSpokes-1);
RadialTraj = (Nbins/2).*(-(-1).^ss).*(rr./Nbins-0.5).*exp(1j.*(pi().*ss./NSpokes));  
RadialTraj = RadialTraj(:);
M = length(RadialTraj);
RadialTraj = unique(RadialTraj); % removing multiple samples of the origin.
GridCoordinates(:,1) = imag(RadialTraj);
GridCoordinates(:,2) = real(RadialTraj);

figure('Name','Sampling trajectory','NumberTitle','off');
plot(GridCoordinates(:,1),GridCoordinates(:,2),'.');
title(['Radial Trajectory with N_{spokes}= ',num2str(NSpokes),', M= ',num2str(M)]);
grid on; grid minor; axis([-0.04 0.02 -0.04 0.02].*Nbins/2)
axis square;

end

