% function [ok,err] = equals(v1,v2,eps)
%
% Compare two multi-dimensional arrays.
%
% v1,v2    multi-dimensional arrays to compare
% eps      Threshold. Default 1.e-12
%
% The comparison is performed element by element. The arrays are considered
% equal if no element in the difference v1-v2 exceeds the threshold.
%
% Returns 1 if the arrays are equal and 0 otherwise.
%
% Yoel Shkolnisky 03/02/03

function [ok,err] = equals(v1,v2,eps)

if nargin<3
    eps = 1.e-12;
end

res = 1;

% Verify that v1 and v2 have the same dimensions
s1 = size(v1);
s2 = size(v2);

if sum(abs(s1-s2))>0
    error('v1 and v2 do not have same dimensions'); %not the same dimensions
end

err=norm(v1(:)-v2(:))/norm(v1(:));
if err > eps
        res=0;
end
ok=res;