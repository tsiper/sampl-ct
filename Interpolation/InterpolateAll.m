function [ ppSinogram ,PSNR_Vals,InterpMethods] = InterpolateAll( Sinogram, N, theta, WindowSize, B, R, W, ppSinogram0 )
%INTERPOLATEALL Interpolates to the Pseudo-Polar grid using all the available
%options

% Defining the standard interpolation options
InterpMethods = {'nearest','linear','pchip','spline','subspace'};

% Preallocating the sinogram vector
ppSinogram = cell(1,length(InterpMethods));

% Interpolating using the different methods
for i=1:length(InterpMethods)-1
    ppSinogram{i}   = InterpolateSinogram( Sinogram, N, theta ,InterpMethods{i});
end

% Interpolating using our method
ppSinogram{i+1} = InterpolateWindow( Sinogram, N, theta, WindowSize, B, R, W );

% Equalizing and measuring psnr values
PSNR_Vals = zeros(length(ppSinogram),1);
for i=1:length(ppSinogram)
   ppSinogram{i} = EqualizeImage(ppSinogram{i},ppSinogram0);
   PSNR_Vals(i)  = psnr(trimcols(ppSinogram{i},N),trimcols(ppSinogram0,N));    
end

end

