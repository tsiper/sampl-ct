function [y] = Diric_m(x,m)
% Dm - Dirichlet Kernel
    q = sin(pi*x/m);
    p = sin(pi*x);
    I = find(q > 1e-14);
    J = 1:length(x);
    J(I) = [];
    y = zeros(length(x));
    y(I) = p(I)./(m.*q(I));
    y(J) = sign(cos(x(J)));
    
end