% Test the inversion of the PPFT using the convolution approach.
%
% Yoel Shkolnisky, December 2010.

n=64;
A=rand(n,n,n);
pp=ppft3(A);
tic;
B=fippft3(pp,1.0e-8,100,1);
err=norm(B(:)-A(:))/norm(B(:));
t=toc;
fprintf('err=%e\n',err);
fprintf('time=%5.2f secconds\n',t);
