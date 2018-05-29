function [ rot_matrix_x ] = x_rotate3D( teta )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

rot_matrix_x = [1 0 0 ; 0 cosd(teta) -sind(teta); 0 sind(teta) cosd(teta)];

end

