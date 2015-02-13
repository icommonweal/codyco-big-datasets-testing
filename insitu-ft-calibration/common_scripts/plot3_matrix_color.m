function [ output_args ] = plot3_matrix_color( input_args, plot_color )
%PLOT3_MATRIX version of plot3 that takes an input a single matrix
if( size(input_args,2) == 3) 
    x = input_args(:,1);
    y = input_args(:,2);
    z = input_args(:,3);
end
if( size(input_args,1) == 3) 
    x = input_args(1,:);
    y = input_args(2,:);
    z = input_args(3,:);
end
output_args = plot3(x,y,z,'Color',plot_color);

end

