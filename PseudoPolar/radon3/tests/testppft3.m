% function testppft3
%
% Test the functions ppft3_ref, ppft3, and radon3.
%
% Yoel Shkolnisky 20/05/2013

function testppft3

eps = 1.e-14;

sz=[4 8 16 20 32 40 64 100 128 200 256];

for n=sz
    % Test the function ppft3 by comparing it to ppft3_ref.
    im = rand(n,n,n);
    tic;
    pp = ppft3_ref(im);
    t1=toc;
    tic;
    pp2 = ppft3(im);
    t2=toc;
    [ok,err]=equals(pp,pp2,eps);
    reportResult('Test ppft3',n,ok,err,t1,t2);
end

% Test the function Radon3 by comparing it to the reference function
% slowRadon3.
im = rand(4,4,4);
tic;
rr = slowradon3(im);
t1=toc;
tic;
rr2 = Radon3(im);
t2=toc;
[ok,err]=equals(rr,rr2,eps);
reportResult('Test radon3',4,ok,err,t1,t2);



function reportResult(testMsg,n,res,err,t1,t2)
if (res)
    str = 'OK';
else
    str = 'FAIL';
end
fprintf('%s \t n=%d \t %s \t err=%e \t t1=%7.5f(secs) \t t2=%7.5f(secs)  speedup=%5.2f\n',testMsg,n,str,err,t1,t2,t1/t2);
