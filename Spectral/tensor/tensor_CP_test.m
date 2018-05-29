mat = [1 2 3; 4 5 6; 7 8 9];
%mat = repmat(mat,50,50);
cell_mat = cell(3,1);

for i = 1:3 
    cell_mat{i,1} = mat;
end

[Denoising_matrix_CP , Denoising_matrix , CP_Rank] = Denoising_CP( cell_mat , 3 , 'SimpleNoise' , 0.03 , -1, 30); 

