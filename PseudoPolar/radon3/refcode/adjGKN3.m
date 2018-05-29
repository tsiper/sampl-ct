% function v = adjGKN3(w,k)
%
% Computes the adjoint of the operator GKN3 for the vector w and the row k.
% For a vector w of length n:
% 		GKN3(w) = U_{3n+1,n+1}(F^{2k/n}_{3n+1}(E_{n,3n+1}(F^{-1}(w)))) 
% Hence,
% 		adjGKN3 = adj F^{-1} \circ adj E_{n,3n+1} \circ adj F^{2k/n}_{3n+1} \circ adj U_{3n+1,n+1}
%
% GKN3 maps a vector of length n into a vector of length n+1. Therefore,
% adjGKN3 maps a vector of length n+1 into a vector of length n. The length
% of the input vector w is n+1 (n even) and the length of the output vector
% is n.
% 
% w    The sequence to resample. Can be of odd or even length.
% k    The row whose transform is computed.
%
% Returns the adjoint of GKN3 for the sequence w and the row k.
% See thesis' final version for more information.
%
% See Also GKN3.
% 
% Yoel Shkolnisky 26/02/03

function v = adjGKN3(w,k)
n = length(w)-1;
if (mod(n,2)==1)
   error('Length of input vector must be n+1 (n even)');
end

v = adjInvF(adjE(adjCfrft(adjU(w),2*k/n)));

%%%%%%%%%%%%%%%     Subroutines   %%%%%%%%%%%%%%%%%%
function z=adjU(v)
% Adjoint of the truncation operator.
% The length of the vector v is assumesed to be n+1 with n even.
n = length(v)-1;
z = zeros(1,3*n+1);
z(n+1:2*n+1) = v;


function z=adjE(v)
% Adjoint of the extension opreator.
% The length of the vector v is assumesed to be 3*n+1 with n even.
n = (length(v)-1)/3;
z = zeros(1,n);
z = v(n+1:2*n);

function z=adjCfrft(x,alpha)
z = cfrft(x,-alpha);

function z=adjInvF(v)
n=length(v);
z = cfft(v)./n;