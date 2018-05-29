function [ p_new ] = ProjToPositive3( U,p_vec )
% U - unitaric matric where it's columns are orthogonal basis vectors
v1 = U(:,1);
v2 = U(:,2);
v3 = U(:,3);
p = reshape(p_vec,3,1);

C12 = cross(v1,v2);
res = C12'*p;
if res<0
    p12_new = p-(C12'*p)/norm(C12)^2*C12;
else
    p12_new = p;
end

C13 = cross(v3,v1);
res = C13'*p12_new;
if res<0
    p13_new = p12_new-(C13'*p12_new)/norm(C13)^2*C13;
else
    p13_new = p12_new;
end

C23 = cross(v2,v3);
res = C23'*p13_new;
if res<0
    p_new = p13_new-(C23'*p13_new)/norm(C23)^2*C23;
else
    p_new = p13_new;
end

end

