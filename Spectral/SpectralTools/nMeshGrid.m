function [ ndots_matrice ] = nMeshGrid( coordinates_matrice )
% calculates meshgrid for size(coordinates_matrice) and returns vectors of
% dim n.
% - coordinates_matrice: MXN matrice, where M is the length of the
% coordinate vectors and N is their dimension;
[M, N] = size(coordinates_matrice);
ndots_matrice = zeros(M^N, N);
for n=1:N
    for m=1:M
        ndots_matrice(1+(m-1)*M^(N-n):m*M^(N-n),n)= repmat(coordinates_matrice(m,n),M^(N-n),1);
    end
    ndots_matrice(:,n)=repmat(ndots_matrice(1:M^(N-n+1),n),M^(n-1),1);
end
ndots_matrice = shiftdim(ndots_matrice,-1);
end

