function [ force_matrix ] = static_force_matrix( mass )
%static_wrench_matrix Return the 3x3 matrix such that maps the gravity 
% to the static wrench measured by the sensor
force_matrix = mass*eye(3);

end

