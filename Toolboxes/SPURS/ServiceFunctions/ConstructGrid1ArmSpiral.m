function [ GridCoordinates ] = ConstructGrid1ArmSpiral( sqrtN, M, Ns)
% This fucnion returns a constatnt speed spiral trajectory.
% The trajectory has a support of a circle with radius=sqrt(N)/2.
% Ns = sqrt(M/pi()) yields a close to uniform density spiral.

A = sqrtN/2;
R = sqrt(((0:M-1)./M));
w = 2*pi()*Ns*R;
GridCoordinates = transpose(A*[R.*cos(w); R.*sin(w)]);

figure('Name','Sampling trajectory','NumberTitle','off');
plot(GridCoordinates(:,1),GridCoordinates(:,2),'.-');
title(['Spiral Trajectory with 1 arms, N= ',num2str(sqrtN^2),', M= ',num2str(M)]);
grid on; grid minor; axis([-0.04 0.02 -0.04 0.02].*sqrtN)
axis square;

end

