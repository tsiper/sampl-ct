function [ col_cm,row_cm ] = find_COM( B )
%Finds coordinates of COM
B(B<10^-7)=0;
B(B>0)=1;
[m,n] = size(B);
M= sum(sum(B));
p_cols= 1:n ;
p_rows= 1:m ;
[col , row ]= meshgrid(p_cols, p_rows); %X is cols, Y is rows

col_cm =round( sum(sum(B.*col)) / M);
row_cm =round( sum(sum(B.*row)) / M);

end

