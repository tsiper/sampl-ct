function [ App ] = BuildPPMtx( n )
%BUILDPPMTX This function builds the PP matrix, for computing the pseudo-polar
%transform, according to equation 2 from Shkolnitzky (I) paper:
%
%              n/2-1  n/2-1
% I_hat(k,l) =  SUM    SUM    I(u,v)*exp(-2j*pi*(-2kl/n*u + k*v))
%             u=-n/2  v=-n/2
%

% We take twice the size of the image as the paramterization of the
% measurements.
m = 2*n+1;

% Initializing the matrix
App = zeros( 2*(n+1)*m , n^2 );

% Looping over all the rows of App, corresponding to the measurements
for i = 1:size(App,1) / 2
    
    % Setting the k,l indices according to the vectorization
    k = mod((i-1),m)   - n;        % Columns of the PP transform
    l = floor((i-1)/m) - n/2;      % Rows of the PP transform
    
    % Looping over all the columns of App, corresponsing to pixels in the image
    for j=1:size(App,2)

        % Setting the i,j indices according to the vectorization
        u = mod((j-1),n)   - (n/2); % Columns of the original image
        v = floor((j-1)/n) - (n/2); % Rows of the original image

        % Finally - computing the corresponding coefficient
        w_u = -2*l/n*k;
        w_v = k;
        % % This is what should have been the formula
        % App(i,j) = exp(-2j*pi/m * (w_u*u+w_v*v) );
        
        App(i,j) = -exp(-2j*pi/m * (w_u*u+w_v*v - k/2) );

        % Translating the rotation to new u and v coordinates
        vt = -u -1;
        ut = v;
        
        % Calculating for the rotated image again with a shift of indices
        App((n+1)*m+i,j) = -exp(-2j*pi/m * (w_u*ut+w_v*vt - k/2 ) );

    end
end


end

