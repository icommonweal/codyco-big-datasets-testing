function [ force_regressor ] = static_force_regressor( )
%static_wrench_regressor Return the 18x1 matrix such that maps 
force_regressor = zeros(18,1);
force_regressor(1,1) = 1;
force_regressor(8,1) = 1;
force_regressor(15,1) = 1;

%check
tol = 1e-5;


end

