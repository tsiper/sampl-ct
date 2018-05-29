function [ x_fix ] = ShiftImage( x,x_ref )
%SHIFTIMAGE Shifts x, so it will register on x_ref

x_fix = circshift(x,FindShift(x,x_ref));

end

