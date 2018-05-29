function [ y_pp ] = PP_Downscale( y_pp_pad )

N = (size(y_pp_pad,1)-1)/4;
DecFactor = 2;

y_pp_temp = trimcols(y_pp_pad,2*N+1);
D = DecOperator(DecFactor,'uniform');
inds = D(ones(1,4*N+2));
y_pp = zeros(2*N+1,2*N+2);

j=1;
for i=1:length(inds)
    if inds(i)
        y_pp(:,j) = y_pp_temp(:,i);
        j = j+1;
    end
end
end

