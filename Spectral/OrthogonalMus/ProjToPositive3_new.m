function [ p_new ] = ProjToPositive3_new( invP,P,p_vec )
% U - unitaric matric where it's columns are orthogonal basis vectors
p_original = (invP*p_vec);
p_original(p_original<0) = 0;
if norm(p_original-(invP*p_vec))~=0
    p_new = P*p_original;
else
    p_new = p_vec;
end
end

