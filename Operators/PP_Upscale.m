function [ y_pp_pad ] = PP_Upscale( y_pp )

DecFactor = 2;

N = (size(y_pp,1)-1)/2;

y_pp_temp = padcols(y_pp,4*N+1);
D = DecOperator(DecFactor,'uniform');
inds = D(ones(1,4*N+2));
y_pp_pad = zeros(4*N+1,4*N+2);

j = 1;
for i=1:length(inds)
    if inds(i)
        y_pp_pad(:,i) = y_pp_temp(:,j);
        j = j+1;
    end
end
end

