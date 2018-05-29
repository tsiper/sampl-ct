function [ I ] = adj_ppfft( Ipp1, Ipp2 )
%IPPFFT 

% Getting the dimensions of the image
[m,~] = size(Ipp1);
n=(m-1)/2;

I1 = zeros(n,m);
I2 = I1;

for k=-n:n
    q = fliplr(Ipp1(k+n+1,:));
    I1(:,k+n+1) = adjG(k,n,q);
end
I1 = invF2(I1);
for k=-n:n
    q = fliplr(Ipp2(k+n+1,:));
    I2(:,k+n+1) = adjG(k,n,q);
end
I2 = invF2(I2);

% Cropping both results
I1c = I1(:,(n/2+1):(end-n/2-1));
I2c = rot90(I2(:,(n/2+1):(end-n/2-1)),1);

I = (I1c+I2c)/n/m;

end

function [w] = adjG(k,n,q)
    q = q(:);
    a = 2*k / n;
    % Adjoint psuedo-polar
    temp = adjFam(q,a);
    w    = F1(temp);
end 