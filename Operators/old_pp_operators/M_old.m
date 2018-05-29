function [ y ] = M_old( y )
%M Preconditioning operator for the conjugate-gradient case
% The column size of y must be odd!

m = size(y,1);

if mod(m,2)
    % When dealing with the odd case
    n = (m-1)/2;
    filt = abs(-n:n)'/m;
    filt(n+1) = 1/(2*m);
else
    filt = abs(-(m-1)/2 : (m-1)/2)'/m;
end

for i=1:size(y,2)
    y(:,i) = y(:,i).*filt;
end

end