function [ wrench_matrix ] = static_wrench_matrix( pi )
%static_wrench_matrix Return the 6x3 matrix such that maps the gravity 
% to the static wrench measured by the sensor
mass = pi(1);
mcog = pi(2:4);
wrench_matrix = zeros(6,3);
wrench_matrix(1:3,:) = mass*eye(3);
wrench_matrix(4:6,:) = skew(mcog);

end

