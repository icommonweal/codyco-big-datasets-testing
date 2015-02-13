function [ output_args ] = scatter3_matrix( input_args )
%PLOT3_MATRIX version of scatter3 that takes an input a single matrix
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
output_args = scatter3(x,y,z,'red');

end
