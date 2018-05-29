% function testAdj
%
% Tests the functions adjPPFT and adjRadon.
% Check if <PPFT(A),B> = <A,adjPPFT(B)> for various A and B. 
% <A,B> stands for the inner-product of A and B.
%
% Legend for results format:
%   n      - matrix size
%   a      - <PPFT(A),B>
%   b      - <A,adjPPFT(B)>
%   error  - absolute error a-b
%   Rel a  - Relative error (a-b)/a
%   Rel b  - Relative error (a-b)/b
%
% See also PPFT, adjPPFT, PsuedoRadon, adjRadon.
%
% Yoel Shkolnisky 9/2/02

function testAdj;

for l=1 %l=1:2
   if (l==1)
      disp('******** START Test of adjPPFT **********');
   else     
      disp('******** START Test of adjRadon **********');
   end;
        
    for k=[4,8,16,32,64]%,128,256,512]
        A = rand(k);
        B1 = rand(2*k+1,k+1);
        B2 = rand(2*k+1,k+1);
      
      if (l==1)
         [r1,r2] = ppft(A);
        a = ip([r1 r2],[B1 B2]);
            b = ip(A,adjPPFT(B1,B2));
      else
         [r1,r2] = Radon(A);
        a = ip([r1 r2],[B1 B2]);
            b = ip(A,adjRadon(B1,B2));
      end;

        disp (strcat('n = ',num2str(k)));
        disp (strcat('a =  ',num2str(a)));
        disp (strcat('b =  ',num2str(b)));
        disp (strcat('Error: ',num2str(a-b)));
        disp (strcat('Rel a =  ',num2str((a-b)/a)));
        disp (strcat('Rel b =  ',num2str((a-b)/b)));
        disp ('*********************');
      disp (' ');
   end;
end;
