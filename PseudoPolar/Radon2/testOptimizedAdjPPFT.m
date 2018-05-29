% function testOptimizedAdjPPFT
%
% Tests the function OptimizedAdjPPFT.
% Check if <optimizedPPFT(A),B> = <A,optimizedAdjPPFT(B)> for various A and B. 
% <A,B> stands for the inner-product of A and B.
%
% Legend for results format:
%   n      - matrix size
%   a      - <optimizedPPFT(A),B>
%   b      - <A,optimizedAdjPPFT(B)>
%   error  - absolute error a-b
%   Rel a  - Relative error (a-b)/a
%   Rel b  - Relative error (a-b)/b
%
% Yoel Shkolnisky 22/10/01

function testOptimizedAdjPPFT

for k=[4,8,16,32,64]%,128,256,512]
    A = rand(k);
    B1 = rand(2*k+1,k+1);
    B2 = rand(2*k+1,k+1);
    
    [r1,r2] = optimizedPPFT(A);
    a = ip([r1 r2],[B1 B2]);
    b = ip(A,OptimizedAdjPPFT(B1,B2));

    disp (strcat('n = ',num2str(k)));
    disp (strcat('a =  ',num2str(a)));
    disp (strcat('b =  ',num2str(b)));
    disp (strcat('Error: ',num2str(a-b)));
    disp (strcat('Rel a =  ',num2str((a-b)/a)));
    disp (strcat('Rel b =  ',num2str((a-b)/b)));
    disp ('*********************');
   disp (' ');
end;
