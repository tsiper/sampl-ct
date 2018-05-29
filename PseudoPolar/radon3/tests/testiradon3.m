% function testIRadon3
%
% Tests the function iRadon3.
%
% Tests the correctness and performance of the inversion algorithm of the 3-D discrete Radon transform.
% The function tests also the inversion of the 3-D pseudo-polar Fourier transform since 
% inverting the 3-D discrete Radon transform requirers inverting the 3-D pseudo-polar Fourier transform.
%
% Yoel Shkolnisky 01/03/03

function testiradon3

% test No.1: All ones cube 8x8x8
% Inversion should converge
a = ones(8,8,8);
pp = radon3(a);
runTest(pp,1.e-2,20,a,'1 - All ones 8x8x8',1,1);

% test No.2: Random cube [0,255] 8x8x8
% Inversion should converge
a = fix(255*(randn(8,8,8)));
pp = radon3(a);
runTest(pp,1.e-2,20,a,'2- Random cube [0,255] 8x8x8',1,1);

% test No.3: Random cube [0,255] 8x8x8
% Inversion should NOT converge (required error too small)
a = fix(255*(randn(8,8,8)));
pp = radon3(a);
runTest(pp,1.e-7,10,a,'3 - Random cube [0,255] 8x8x8 (error too small)',1,1);

% test No.4: Random cube [0,255] 8x8x8
% Inversion should NOT converge (not enough iterations)
a = fix(255*(randn(8,8,8)));
pp = radon3(a);
runTest(pp,1.e-4,5,a,'4 - Random cube [0,255] 8x8x8 (not enough iterations)',1,1);


%test No.5: Random cube 8x8x8 with real values from [0,1]
a = randn(8,8,8);
pp = radon3(a);
runTest(pp,1.e-2,100,a,'5 - Real random cube 8x8x8 from [0,1]',1,0);

%test No.8: Random cube [0,255] 8x8x8 (run quietly)
a = fix(255*rand(8,8,8));
pp = radon3(a);
runTest(pp,1.e-2,10,a,'6 - Random cube 8x8x8 [0,255] (quiet)',0,1);


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
% trueimage    True if the input represents the discrete Radon transform of an integer image.
%              In this case it is possible to exactly inverting the transform
%              by truncating any floating parts of the result.
function runTest(pp,ErrTol,MaxIts,ref,description,verbose,trueimage)

fprintf('Test name : %s\n',description);
fprintf('Requested error = %-2.5e\n',ErrTol);
fprintf('Max number of iterations = %d\n\n',MaxIts);
if verbose
    fprintf('Inversion log:\n');
   fprintf('--------------\n');
end
tic;
[Y,flag,res,iter] = iradon3(pp,ErrTol,MaxIts,verbose);
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
diff = abs(Y-ref);
fprintf('Absolute error = %e\n',max(diff(:)));

if trueimage
    diff = abs(round(Y)-ref);
    fprintf('Absolute error of reconstructed image = %-1.3f\n',max(diff(:)));
end
   
fprintf('Computation time = %-3.2f secs\n',t);
fprintf('--------------------------------------------------------\n\n\n');
