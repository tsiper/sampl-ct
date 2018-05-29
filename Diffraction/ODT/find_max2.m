function [ max_ind , max_val ] = find_max( vec, cols )
abs_vec= abs(vec);
max_ind = cols/2+5;
for k = cols/2+6:cols %not including DC
    if abs_vec(k)> abs_vec(max_ind)
        max_ind= k;
    end
end

max_val= abs_vec(max_ind);

end

