function [ r,s ] = Proj_P( p,q )
%PP Projects the x and y derivatives of p and q to r and s according to Amir
%Beck's paper

% The dimensions of the objective
m = size(q,1);
n = size(p,2);

% padding the dimensions of p q to (m x n)
p = [zeros(1,n);p];
q = [zeros(m,1),q];

% Calculating the normalizing matrix
W = max(1,sqrt(p.^2+q.^2));

% Normalizing
r = p./W;
s = q./W;

% Removing the unnecessary columns
r = r(2:end,:);
s = s(:,2:end);

end

