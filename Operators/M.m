function [ y ] = M( y )
%M Preconditioning operator for the psuedo-polar radon transform

% Extracting the dimensions of the PP sinogram
[M,N] = size(y);
n = (M-1)/2;

% If y is real, than we want to work on the Fourier transform of the columns. Otherwise, do nothing
RealFlag = 0;
if isreal(y)
	RealFlag = 1;
	y = F1(y);
end

% Bulding the Ramp filter for odd and even cases
if mod(M,2)
    % When dealing with the odd case
    filt = abs(-n:n)'/(2*n);
    filt(n+1) = 1/(5.46*2*n);
else
	warning('The Pseudo-Polar sinogram has an even number of elements in its columns');
    filt = abs(-n:n)'/(2*n+1);
end

% Applying the filter, while nullifying the redundant columns 
% (twice 45 degrees, i=N/2, and -45,135 redundancy, i=N)
for i=1:N
    y(:,i) = y(:,i).*filt;
%     if ((i==N/2)||(i==N))&&(i>1)
%         y(:,i) = zeros(M,1);
%     end
end

% Returning to a real number if RealFlag is 'true'
if RealFlag
	y = invF1(y);
end

end