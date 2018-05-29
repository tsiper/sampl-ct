function [ KR ] = KhatriRao( A , B )

N = size(A,2);
if N~=size(B,2)
    error('ColumnMismatch','Input matrices must have the same number of columns');
end

KR = arrayfun(@(ii) kron(A(:,ii),B(:,ii)),1:N,'UniformOutput',0);
KR = cell2mat(KR);


end

