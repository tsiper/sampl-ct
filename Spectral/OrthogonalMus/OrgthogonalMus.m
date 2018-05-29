function [ ort_mus,invP_operator, P_operator, eigs ,U] = OrgthogonalMus(Mu )
% based on calculation described in "An algorithm for constrained one-step inversion of spectral CT data
%Barber R Sidky E Schmidt T Pan X", Pg. 18-19.
%   Mu   -                          matrix where its rows are the atenuation coeff. of different
%                                   materials.
%   ort_mus     -                   matrix where its rows are the orthogonalized atenuation coeff.
%   invP_operator, P_operator  -    as defined in the paper.

% according to the paper: mu_effective = SUM_m (R*(mu*invP)*(P*w))=SUM_m (R*mu_hat*w_hat) 
M = Mu*Mu';
[U,D] = eig(M);
eigs = diag(D);
P_operator = sqrt(D)*U';
invP_operator = U*1/sqrt(D);
ort_mus = Mu'*invP_operator;

end

