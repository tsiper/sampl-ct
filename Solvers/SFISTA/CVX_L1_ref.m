function [ x_out ] = CVX_L1_ref( b, A, Lambda, D, Dt, Type )
%CVX_L1 - Solves L1 regularization problem

[M, N] = size(A);

switch lower(Type)
    case 'synthesis'
        cvx_begin
            variable z(N,1)
            minimize square_pos(norm(A*D(z) - b)) + Lambda*norm(z,1)
        cvx_end
        
        x_out = D(z);
    case 'analysis'
        cvx_begin
            variable x(N,1)
            minimize square_pos(norm(A*x - b)) + Lambda*norm(Dt(x),1)
        cvx_end
        
        x_out = x;
end

