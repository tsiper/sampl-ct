function n=verifyPP(pp)
%
% Verify that the input image is of size 3x(3n+1)x(n+1)x(n+1) with n even.
% If so, the function return n. Otherwise, the function terminates with
% proper error message.
%
% This function was factored from the functions in ../ppft
%
% Yoel Shkolnisky, December 2010.

ERR = 'pp must be of size 3x(3n+1)x(n+1)x(n+1)';
s = size(pp);

if length(s)~=4
    % The input array is not a 4-D array
    error(ERR);
end

if s(1)~=3
    % Length of first dimension is not 3
    error(ERR);
end

if (s(3)-s(4))~=0
    % The last two dimensions of pp should be n+1
    error(ERR);
end

n = s(3)-1;
if (mod(n,2)~=0)
    % n is not even
    error(ERR);
end

if s(2)~=3*n+1
    % Second dimension is not 3n+1
    error(ERR);
end
