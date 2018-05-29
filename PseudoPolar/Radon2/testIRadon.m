% function testIRadon
%
% Tests the functions iRadon and ippft.
%
% Tests the correctness and performance of the inversion algorithm of the discrete Radon transform.
% The function tests also the inversion of the pseudo-polar Fourier transform since 
% inverting the discrete Radon transform requirers inverting the pseudo-polar Fourier transform.
%
% Yoel Shkolnisky 21/12/02

function testIRadon

% test No.1: magic square of size 64x64
% Inversion should converge
a = magic(64);
[pp1,pp2] = Radon(a);
runTest(pp1,pp2,1.e-2,20,a,'1 - Magic square 64x64',1,1);

% test No.2: magic square of size 64x64
% Inversion should NOT converge (required error too small)
a = magic(64);
[pp1,pp2] = Radon(a);
runTest(pp1,pp2,1.e-7,10,a,'2 - Magic square 64x64 (error too small)',1,1);

% test No.3: magic square of size 64x64
% Inversion should NOT converge (not enough iterations)
a = magic(64);
[pp1,pp2] = Radon(a);
runTest(pp1,pp2,1.e-4,5,a,'3 - Magic square 64x64 (not enough iterations)',1,1);


%test No.4: magic square of size 128x128
a = magic(128);
[pp1,pp2] = Radon(a);
runTest(pp1,pp2,1.e-1,10,a,'4 - Magic square 128x128',1,1);

%test No.5: random 64x64 image with real values from [0,1]
a = rand(64);
[pp1,pp2] = Radon(a);
runTest(pp1,pp2,1.e-2,10,a,'5 - Real random matrix 64x64 from [0,1]',1,0);

%test No.6: random 64x64 image with integer values from [0,255]
a = fix(rand(64)*256);
[pp1,pp2] = Radon(a);
runTest(pp1,pp2,1.e-2,10,a,'6 - Integer random matrix 64x64 from [0,255]',1,1);

%test No.7: Lenna 256
a = double(imread('C:\Private\University\Images\lena256.bmp'));
[pp1,pp2] = Radon(a);
runTest(pp1,pp2,1.e-2,10,a,'7 - Lena 256',1,1);

%test No.8: magic square 32x32. Run quietly
a = magic(32);
[pp1,pp2] = Radon(a);
runTest(pp1,pp2,1.e-2,10,a,'8 - Magic square 32x32 (quiet)',0,1);


%%%%%%%%%%%%%%%
% Sub functions
%%%%%%%%%%%%%%%

% Execute a single inversion test.
% pp1,pp2      The Radon sectors.
% ErrTol       Residual error required from the inversion algorithm.
% MaxIts       Number of iterations of the inversion algorithm.
% ref          The original image. Used as a reference to check the absolute error.
% description  Test description for printing purposes.
% verbose      If true, print the inversion log.
% trueimage    True if the input represents the discrete Radon transform of an integer image.
%              In this case it is possible to exactly inverting the transform
%              by truncating any floating parts of the result.
function runTest(pp1,pp2,ErrTol,MaxIts,ref,description,verbose,trueimage)

fprintf('Test name : %s\n',description);
fprintf('Input size = %d\n',(length(pp1)-1)/2);
fprintf('Requested error = %-2.5e\n',ErrTol);
fprintf('Max number of iterations = %d\n\n',MaxIts);
if verbose
    fprintf('Inversion log:\n');
   fprintf('--------------\n');
end
tic;
[Y,flag,res,iter] = iRadon(pp1,pp2,ErrTol,MaxIts,verbose);
t = toc;
if ~flag
   str = 'CONVERGED';
else
   str = 'DID NOT CONVERGE';
end

fprintf('\nResults:\n');
fprintf('---------\n');
fprintf('Inversion %s\n',str);
fprintf('Residual error %-2.5e at iteration no. %d\n',res,iter);
fprintf('Absolute error = %-2.5f\n',max(max(abs(Y-ref))));

if trueimage
   fprintf('Absolute error of reconstructed image = %-1.3f\n',max(max(abs(round(Y)-ref))));
end
   
fprintf('Computation time = %-3.2f secs\n',t);
fprintf('--------------------------------------------------------\n\n\n');
