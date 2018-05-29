function [ mu,b ] = PowerMethod(A,InputSize)
%calculates max eigenvalue of A. A has to be a square operator
% InputSize is an array that contains the dimensions for the input size for the
% functions A()

b = randn(InputSize);
b = b./norm(b,'fro');

mu=0;
err = 1e-3;
normAb = norm(A(b),'fro');
diff = norm(A(b)-mu*b,'fro')/normAb;

wb = MyWaitbar(0,'Computing Lipschitz Constant - Power Method...');
iter = 0;
max_iters = 50;

while (diff>err) && (iter < max_iters)
    % Updating b
    Ab = A(b);
    b = Ab./norm(Ab(:));
    
    
    % Updating mu for this iteration
    mu = b(:)'*Ab(:) / norm(b(:))^2;
%     L = b'*A'*(A*b);
    diff = norm(A(b)-mu*b,'fro')/normAb;
    
    iter = iter+1;
    MyWaitbar(iter/max_iters,wb,sprintf('Power Method, Iteration %d, Diff=%.4s',iter,diff));
    
end

close(wb);

end