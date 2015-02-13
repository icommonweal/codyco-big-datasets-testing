function [ wrench_regressor ] = static_wrench_regressor( )
%static_wrench_regressor Return the 18x4 matrix such that maps the gravity 
% to the static wrench measured by the sensor

wrench_regressor = zeros(18,4);
wrench_regressor(1,1) = 1;
wrench_regressor(8,1) = 1;
wrench_regressor(15,1) = 1;
wrench_regressor(12,2) = 1;
wrench_regressor(17,2) = -1;
wrench_regressor(16,3) = 1;
wrench_regressor(6,3) = -1;
wrench_regressor(5,4) = 1;
wrench_regressor(10,4) = -1;

%check
tol = 1e-5;


end

