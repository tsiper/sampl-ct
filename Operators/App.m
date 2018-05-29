function [ R ] = App( I , Mproj )
%A Applies the forward pp-radon transform 

% Check that the image is square
if size(I,1) ~= size(I,2)
    error('The input image is not square');
end
n = size(I,1);

% If not m is selected take the default one
if nargin == 2
    PadSize = floor((Mproj-n)/2);
    if PadSize > 0
        I = padarray(I,[PadSize,PadSize]);
    else
        PadSize = 0;
    end
else
    PadSize = 0;
end

% Performing the PP transform on each dimension of the image
% % Ipp1 = ppfft(I);
% % Ipp2 = ppfft(rot90(I,3));

[Ipp1, Ipp2] = OptimizedPPFT(I);
% [Ipp1, Ipp2] = ParOptimizedPPFT(I);


% R  = [invF1(Ipp1),invF1(Ipp2)];
% R  = real(R);
R = invF1([fliplr(Ipp2), Ipp1]);

% Trimming back according to our padding
R = F1( R( 2*PadSize+1 : 2*PadSize+(2*n+1) , : ) );

end

