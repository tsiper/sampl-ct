function [W] = HammingWindow(t,theta,WindowSize)
% Implements a radial Hamming window

[THETA,T] = meshgrid(theta,t);
% [T,THETA] = meshgrid(t,theta);
R = sqrt(T.^2+(THETA/pi).^2);

W = zeros(size(R));
W(R<WindowSize/2) = 0.54 + 0.46*cos(2*pi*R(R<WindowSize/2) / WindowSize) ;

end
    
