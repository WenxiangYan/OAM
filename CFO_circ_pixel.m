function [ mask ] = CFO_circ_pixel( N_row, N_col, offset_row, offset_col, radius )
%CFO_circ_pixel Draw a solid circle
[X,Y]=meshgrid(1:N_col,1:N_row);
center_row=(N_row-1)/2+1+offset_row;
center_col=(N_col-1)/2+1+offset_col;
mask=sqrt((X-center_col).^2+(Y-center_row).^2)<radius;

end

