function [ I ] = App_T( R )
%ADJA Summary of this function goes here
%   Detailed explanation goes here
% 
% R1 = R(:,1:end/2);
% R2 = R(:,end/2+1:end);

% Ipp1 = F1(R1);
% Ipp2 = F1(R2);
% 

Mproj = size(R,2);
n     = (size(R,1)-1)/2;
n_pad = Mproj -1;

% Padding the array back to its size if necessary
if Mproj > n+1
    % The same pad size from App
    PadSize = floor((n_pad-2*n-1)/2);
    R = invF1( padarray( F1(R) , PadSize ) );
else
    PadSize = 0;
end

pp1 = R(:,1:end/2);
pp2 = R(:,end/2+1:end);

% I = real(precondAdjPPFT(pp1,pp2)*2);
I = OptimizedadjPPFT(pp2,fliplr(pp1));


% And now trimming back to the original size
I = I( PadSize/2+1 : PadSize/2+n , PadSize/2+1 : PadSize/2+n);

% Normalizing back
[m,~] = size(pp1);
n=(m-1)/2;
I = I / (n+1)^2;

end

