function [ Res ] = ppfft( I )
%PPFFT Summary of this function goes here
%   Detailed explanation goes here

[n1,n2] = size(I);
if (n1 ~= n2) || (mod(n1,2) ~= 0)
    error('Please insert a square image I~(NxN), where N is even');
end
n = n1;
% m = 2*n+1;
m = 2*n+1;

% We calculate pfft for the first dimension only
Id_hat = F2( E_pad(I,m));
% Id_hat = [Id_hat,Id_hat(:,1)];
% Id_hat = [Id_hat;Id_hat(1,:)];
    
% Res = zeros(m,n+1);
Res = zeros(m,n+1);

% Step 2 - Calculate F2(E(I))
for k=-n:n
    q   = Id_hat(:,k+n+1);
    w_k = G(k,n,q(:));
    Res(k+n+1,:) = flipud(w_k(:));
end

end

function w = G(k,n,q)
    q = q(:);
    a = 2*k / n;
    temp = invF1(q);
    w    = Fam(temp,a);
end
