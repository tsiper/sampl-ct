function [ detector_pattern ] = pattern_algorithm(factor , phantom )
% this function will find the optimal detector pattern based on a given
% phantom

%phantom = LoadPhantom(128,'brain');
[m , n] = size(phantom);
detector_pattern = zeros(m,n);
phantom_padded = inf + zeros(m+2,n+2);
phantom_padded(2:m+1 , 2:n+1) = phantom;
start_point = [ceil(2+(m-1).*rand(1)) , ceil(2+(n-1).*rand(1))];
number_of_patterns = n * m / factor;
for pattern_num = 1:number_of_patterns
    while (phantom_padded(start_point(1) , start_point(2)) == inf)
         start_point = [ceil(2+(m-1).*rand(1)) , ceil(2+(n-1).*rand(1))];
    end
    for index = 1:factor
            gradient_decent = [abs(phantom_padded(start_point(1),start_point(2)) - phantom_padded(start_point(1)+1,start_point(2))), ...
                           abs(phantom_padded(start_point(1),start_point(2)) - phantom_padded(start_point(1),start_point(2)+1)), ...
                           abs(phantom_padded(start_point(1),start_point(2)) - phantom_padded(start_point(1)-1,start_point(2))), ...
                           abs(phantom_padded(start_point(1),start_point(2)) - phantom_padded(start_point(1),start_point(2)-1)), ...
                           abs(phantom_padded(start_point(1),start_point(2)) - phantom_padded(start_point(1)-1,start_point(2)+1)), ...
                           abs(phantom_padded(start_point(1),start_point(2)) - phantom_padded(start_point(1)-1,start_point(2)-1)), ...
                           abs(phantom_padded(start_point(1),start_point(2)) - phantom_padded(start_point(1)+1,start_point(2)-1)), ...
                           abs(phantom_padded(start_point(1),start_point(2)) - phantom_padded(start_point(1)+1,start_point(2)+1))];
         phantom_padded(start_point(1),start_point(2)) = inf;   
         detector_pattern(start_point(1)-1,start_point(2)-1) = pattern_num;
         expression=find(gradient_decent == min(gradient_decent));
         if(min(gradient_decent)~=inf &&~isnan(min(gradient_decent)))
             switch expression(1)
                   case 1
                        start_point = [start_point(1)+1,start_point(2)];
                   case 2
                        start_point = [start_point(1),start_point(2)+1];
                   case 3
                        start_point = [start_point(1)-1,start_point(2)];   
                   case 4
                        start_point = [start_point(1),start_point(2)-1]; 
                   case 5
                        start_point = [start_point(1)-1,start_point(2)+1];
                   case 6
                        start_point = [start_point(1)-1,start_point(2)-1];
                   case 7
                        start_point = [start_point(1)+1,start_point(2)-1];
                   case 8
                        start_point = [start_point(1)+1,start_point(2)+1];
             end
         end
    end
end
figure()
imagesc(detector_pattern)
title('detectors pattern');
end

