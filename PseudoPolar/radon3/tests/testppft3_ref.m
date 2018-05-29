% function testRadon3
%
% Test the functions radon3, ppft3_ref, and ppft3
%
% Yoel Shkolnisky 20/05/2013

function testppft3_ref

eps = 1.e-14;

% Test the reference function slowppft3.
% According to the Fourier slice theorem for the 3-D discrete radon
% transform, the 1-D inverse FFT of the pseudo-polar Fourier transform
% (slowppft3) along the parameter k should be equal to the discrete 3-D 
% radon transform.
im = rand(2,2,2);
rr = slowradon3(im);
pp = slowppft3(im);
rr2 = icfftd(pp,2);
[ok,err]=equals(rr,rr2,eps);
reportResult('Test slowppft3',ok,err);

% Test the function ppft3_ref by comparing it to slowppft3 (we assume the
% slowppft3 was verified in the pervious step)
%im = rand(4,4,4);
im = rand(2,2,2);
pp = slowppft3(im);
pp2 = ppft3_ref(im);
[ok,err]=equals(pp,pp2,eps);
reportResult('Test 1 ppft3_ref',ok,err);

im = rand(4,4,4);
pp = slowppft3(im);
pp2 = ppft3_ref(im);
[ok,err]=equals(pp,pp2,eps);
reportResult('Test 2 ppft3_ref',ok,err);


function reportResult(testMsg,res,err)
if (res)
    str = 'OK';
else
    str = 'FAIL';
end
fprintf('%s \t %s \t err=%e\n',testMsg,str,err);
