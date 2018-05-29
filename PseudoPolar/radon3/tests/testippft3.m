% function testippft3
%
% Tests the functions ippft3_ref and ippft3.
%
% Tests the correctness and performance of the inversion algorithm of the 3-D pseudo-polar Fourier transform.
%
% Yoel Shkolnisky 20/05/2013

function testippft3

%sz=[4 8 16 20 32 40 64 100 128 200 256];
sz=[4 8 16 20 32 40 64];

fprintf('n \t err_ref \t err_ppft3 \t err_conv \t t_ref \t\t t_ppft3 \t t_conv \t t_ref/t_ppft3 \t t_ref/t_conv\n');
for n=sz
    % Test the function ppft3 by comparing it to ppft3_ref.
    im = rand(n,n,n);
    pp = ppft3(im);
    runTest(n,pp,1.e-8,100,im,0);
end


%%%%%%%%%%%%%%%
% Sub functions
%%%%%%%%%%%%%%%

% Execute a single inversion test.
% pp           The 3-D Radon sectors.
% ErrTol       Residual error required from the inversion algorithm.
% MaxIts       Number of iterations of the inversion algorithm.
% ref          The original image. Used as a reference to check the absolute error.
% description  Test description for printing purposes.
% verbose      If true, print the inversion log.
function runTest(n,pp,ErrTol,MaxIts,ref,verbose)

tic;
[Y1,~,~,~] = ippft3_ref(pp,ErrTol,MaxIts,verbose);
t1=toc;

tic;
[Y2,~,~,~] = ippft3(pp,ErrTol,MaxIts,verbose);
t2=toc;

tic;
[Y3,~,~,~] = fippft3(pp,ErrTol,MaxIts,verbose);
t3=toc;

err_ref=norm(Y1(:)-ref(:))/norm(ref(:));
err_ppft3=norm(Y2(:)-ref(:))/norm(ref(:));
err_conv=norm(Y3(:)-ref(:))/norm(ref(:));

fprintf('%d \t %e \t %e \t %e \t %6.3f \t %6.3f \t %6.3f \t\t %5.2f \t\t %5.2f \n',...
    n,err_ref,err_ppft3,err_conv,t1,t2,t3,t1/t2,t1/t3);

   