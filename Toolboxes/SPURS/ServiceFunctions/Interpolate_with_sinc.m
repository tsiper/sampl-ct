function [ NonUniform_sample_values ] = Interpolate_with_sinc(N, Uniform_grid_points, Uniform_sample_values, NonUniform_grid_points )

n = sqrt(size(Uniform_grid_points,1));
I = (n/N);

[pku_c, pknu_c]=meshgrid(Uniform_grid_points(1:n,2), NonUniform_grid_points(:,2));
[pku_r, pknu_r]=meshgrid(Uniform_grid_points(1:n:end,1), NonUniform_grid_points(:,1));

col_mat = sinc(pku_c-pknu_c);
row_mat = sinc(pku_r-pknu_r);

clear pku_c; clear pknu_c;
clear pku_r; clear pknu_r;

col_mat = col_mat./I./N;
row_mat = row_mat./I./N;

img = reshape(Uniform_sample_values,n,n);

col_mat = col_mat*img;
row_mat = row_mat.*col_mat;
NonUniform_sample_values = sum(row_mat,2);

end

