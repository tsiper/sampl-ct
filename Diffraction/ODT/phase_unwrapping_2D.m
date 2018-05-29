function [ unwrapped_phase ] = phase_unwrapping_2D( wrapped_phase )
%discrete cosine transform (DCT) based un-weighted least squares (UWLS)
%algorithm-
%based on the DCT method from page 200 in
%two dimensional phaseunwrapping/ D.C. Ghiglia and M.D. Pritt

[rows, cols] = size(wrapped_phase);

temp= zeros(rows, cols);
dif_x= zeros(rows, cols);
dif_y= zeros(rows, cols);

%calculate wrapped differnces
%dif_X
for m=1:rows %for the last column we have zeros due to mirror reflection causing replciation
    for n=1:cols-1
         temp(m,n)= wrapped_phase(m, n+1)- wrapped_phase(m, n); 
         dif_x(m,n)= atan2(sin(temp(m,n)), cos(temp(m,n))); %wrap
    end
end

%dif_Y
for m=1:rows-1 %for the last row we have zeros due to mirror reflection causing replciation
    for n=1:cols
         temp(m,n)= wrapped_phase(m+1, n)- wrapped_phase(m, n); 
         dif_y(m,n)= atan2(sin(temp(m,n)), cos(temp(m,n))); %wrap
    end
end

ro_x= zeros(rows, cols);
ro_y= zeros(rows, cols);

%ro_x
for m=1:rows 
         ro_x(m,1)= dif_x(m, 1); %dif_x(m, 0) is zero due to periodicity and mirroring
    for n=2:cols
         ro_x(m,n)= dif_x(m, n)- dif_x(m, n-1); 
    end
end

%ro_y
for n=1:cols
         ro_y(1,n)= dif_y(1, n); %dif_x(0, n) is zero due to periodicity and mirroring
end

for m=2:rows 
    for n=1:cols
         ro_y(m,n)= dif_y(m, n)- dif_y(m-1, n); 
    end
end

%ro
ro= ro_x+ro_y;

result=dct2(ro);

  	for m=1:rows
        for n=1:cols
            if (n==1&&m==1)
                result(m,n)= 0; 
            else
            result(m,n)= result(m,n)/( (2*cos(pi*(m-1)/rows) + 2*cos(pi*(n-1)/ cols) -4));
        end
        end
    end        

unwrapped_phase=idct2(result);
 %%