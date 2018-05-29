% function testRadon
%
% Tests the functions PPFT,Optimized PPFT, and Radon.
% See also PPFT, OptimizedPPPT, Radon.
%
%  Yoel Shkolnisky 9/2/02

function testRadon

eps = 0.00000001;
im = magic(8);

[spp1,spp2] = slowPPFT(im);
[pp1,pp2]   = PPFT(im);
[opp1,opp2] = OptimizedPPFT(im);

if isempty(find(spp1-pp1>eps)) | isempty(find(spp2-pp2>eps))
    disp('PPFT OK');
else    disp('PPFT NOT OK');
end

if isempty(find(spp1-opp1>eps)) | isempty(find(spp2-opp2>eps))
    disp('OptmizedPPFT OK');
else    disp('OptmizedPPFT NOT OK');
end


[sr1,sr2] = slowRadon(im);
[r1,r2]   = Radon(im);

if isempty(find(sr1-r1>eps)) | isempty(find(sr2-r2>eps))
    disp('Radon OK');
else    disp('Radon NOT OK');
end

 
