% Test the inversion of the PPFT using the GPU acceleration
%
% Yoel Shkolnisky, 21/05/2013.

n=64;
A=rand(n,n,n);
pp=ppft3(A);
tic;
B=fippft3_gpu(pp,1.0e-8,100,1);
err=norm(B(:)-A(:))/norm(B(:));
t=toc;
fprintf('err=%e\n',err);
fprintf('time=%5.2f secconds\n',t);
