% function v = adjGKN(w,k)
%
% Computes the adjoint of the operator GKN for the vector w and the row k.
% For a vector w of length n:
% 		GKN(w) = U_{2n+1,n+1}(F^{2k/n}_{2n+1}(E_{n,2n+1}(F^{-1}(w)))) 
% Hence,
% 		adjGKN = adj F^{-1} \circ adj E_{n,2n+1} \circ adj F^{2k/n}_{2n+1} \circ adj U_{2n+1,n+1}
% (where n is the length of the vector w and not of the input to adjGKN)
%
% w    The sequence to resample. Can be of odd or even length.
% k    The row whose transform is computed.
%
% Returns the adjoint of GKN for the sequence x and the row k.
% See thesis' final version for more information.
%
% See Also GKN.
% 
% Yoel Shkolnisky 22/10/01

function v = adjGKN(w,k)
n = length(w)-1;
v = adjInvF(adjE(adjCfrft(adjU(w),2*k/n)));

%%%%%%%%%%%%%%%     Subroutines   %%%%%%%%%%%%%%%%%%
function z=adjU(v)
n = length(v)-1;
z = zeros(1,2*n+1);
z(n/2+1:3*n/2+1) = v;


function z=adjE(v)
n = (length(v)-1)/2;
z = zeros(1,n);
z = v(n/2+1:3*n/2);

function z=adjCfrft(x,alpha)
z = cfrft(x,-alpha);

function z=adjInvF(v)
n=length(v);
z = cfft(v)./n;