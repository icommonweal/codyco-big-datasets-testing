clear

addpath('../quadfit/')
addpath('../common_scripts')
addpath('../plot2svg')

% Stupid workaround for issue http://stackoverflow.com/questions/19268293/matlab-error-cannot-open-with-static-tls
ones(10,10)*ones(10,10);

%input dataset names
datasets_names = {'dataset0.csv' ...
                  'dataset1.csv' ...
                  'dataset2.csv' ...
                  'dataset3.csv' ...
                  'dataset4.csv'};
              
n_datasets = length(datasets_names);

              
drop_at_start = zeros(1,n_datasets);
              
drop_at_end = zeros(1,n_datasets);
              

added_mass = [ 0.0 0.49156 1.0357 1.5797 2.1237 ];
           
       
added_center_of_mass_wrt_leg{1} = [0.0 0.0 0.0];
added_center_of_mass_wrt_leg{2} = [0.0 0.0 0.0];
added_center_of_mass_wrt_leg{3} = added_center_of_mass_wrt_leg{1};
added_center_of_mass_wrt_leg{4} = added_center_of_mass_wrt_leg{1};
added_center_of_mass_wrt_leg{5} = added_center_of_mass_wrt_leg{1};

sampling_time = 0.01;
%% Loading data

%start building structure array
datasets = {};
n_datasets = length(datasets_names);
for i = 1:n_datasets
    datasets{i} = struct();
    datasets{i}.name = datasets_names{i};
    datasets{i}.drop_at_start = drop_at_start(i);
    datasets{i}.drop_at_end = drop_at_end(i);
end

for i = 1:n_datasets
    fprintf(['Loading and synchronizing dataset ' datasets{i}.name '\n']);
    buf = dlmread(datasets_names{i});
    
    datasets{i}.ft_raw = buf(:,1:6);
    datasets{i}.acc_raw = buf(:,7:9);
    %rot_matrix = [-1,0,0;0,-1,0;0,0,1]
    rot_matrix = eye(3);
    
    datasets{i}.acc = (rot_matrix*datasets{i}.acc_raw')';
    datasets{i}.acc_raw = datasets{i}.acc;
end

%% Clean data


%% Subspace estimation
square_side = ceil(sqrt(n_datasets));
figure
for i = 1:n_datasets
    datasets{i}.ft_raw_no_mean = datasets{i}.ft_raw-ones(size(datasets{i}.ft_raw,1),1)*mean(datasets{i}.ft_raw);
    datasets{i}.ft_raw_mean = mean(datasets{i}.ft_raw);
    
    subplot(square_side,square_side,i);
    [U_raw,S_ft_raw,V_raw] = svd(datasets{i}.ft_raw_no_mean,'econ');
    bar((S_ft_raw));
    datasets{i}.ft_raw_projector = V_raw(:,1:3)';
    datasets{i}.ft_raw_projected = (V_raw(:,1:3)'*datasets{i}.ft_raw_no_mean')';
    title('Raw sensor data singular values')

end

subsampling = 1;
% normalize = @(x) x./norm_factor;
% denormalize = @(x) norm_factor.*x;

normalize = @(x) (x-ones(size(x,1),1)*mean(x))./(ones(size(x,1),1)*std(x));
normalize_isotropically = @(x) (x-ones(size(x,1),1)*mean(x))/mean(std(x));
% denormalize = @(x,y,z) (x.*(ones(size(x,1),1)*z)+ones(size(x,1),1)*y);

% normalize data
for i=1:n_datasets
    datasets{i}.ft_raw_projected_norm = normalize(datasets{i}.ft_raw_projected);
end


%% 
%Plotting ellipsoid fitted in raw space

for i=1:n_datasets
    fprintf(['Fitting ft ellipsoid for dataset ' datasets{i}.name '\n']);
    [datasets{i}.p_ft_norm,datasets{i}.ft_proj_norm_refitted]   = ellipsoidfit_smart(datasets{i}.ft_raw_projected_norm(1:subsampling:end,:),datasets{i}.acc(1:subsampling:end,:));
    fprintf(['Fitting acc ellipsoid for dataset ' datasets{i}.name '\n']);
    datasets{i}.p_acc  = ellipsoidfit(datasets{i}.acc(1:subsampling:end,1),datasets{i}.acc(1:subsampling:end,2),datasets{i}.acc(1:subsampling:end,3));
end

figure
for i = 1:n_datasets
    subplot(square_side,square_side,i);
    plot_ellipsoid_im(datasets{i}.p_ft_norm);
    plot3_matrix(datasets{i}.ft_raw_projected_norm(1:subsampling:end,:))
    axis equal
    title('Ellipsoid fitted in FT raw space');
end

% Offset estimation
for i = 1:n_datasets
    [centers,ax] = ellipsoid_im2ex(datasets{i}.p_ft_norm);
    datasets{i}.center_ft_proj = denormalize2(centers',mean(datasets{i}.ft_raw_projected),std(datasets{i}.ft_raw_projected));
    datasets{i}.offset_ft = ((datasets{i}.ft_raw_projector')*datasets{i}.center_ft_proj')'+datasets{i}.ft_raw_mean;
    datasets{i}.ft_raw_no_offset = datasets{i}.ft_raw - ones(size(datasets{i}.ft_raw,1),1)*datasets{i}.offset_ft;
end





%% 
%Plotting ellipsoid fitted in force space

for i=1:n_datasets
    fprintf(['Fitting ft ellipsoid for dataset ' datasets{i}.name '\n']);
    [datasets{i}.p_force,datasets{i}.force_refitted]   = ellipsoidfit_smart(datasets{i}.ft_raw_no_offset(1:subsampling:end,1:3),datasets{i}.acc(1:subsampling:end,:));
end

figure
for i = 1:n_datasets
    subplot(square_side,square_side,i);
    plot_ellipsoid_im(datasets{i}.p_force);
    plot3_matrix(datasets{i}.ft_raw_no_offset(1:subsampling:end,1:3))
    axis equal
    title('Ellipsoid fitted in force space');
end

%% 
% Plotting accelerometer ellipsoid

figure
for i = 1:n_datasets
    subplot(square_side,square_side,i);
    plot_ellipsoid_im(datasets{i}.p_acc);
    plot3_matrix(datasets{i}.acc(1:subsampling:end,:))
    axis equal
    title('Ellipsoid fitter for accelerometer measure');
end

figure
for i = 1:n_datasets
    datasets{i}.mass_estimated = {};
    datasets{i}.mass_estimated_old = {};
end

% 
%% Calibrate force matrix


for j=2:n_datasets

calibration_datasets = 1:j;
validation_datsets = 3:5;

acc_ft = zeros(18+1,1);
cov_ft = zeros(18+1,18+1);

for cal_dat = calibration_datasets 
    for smpl = 1:size(datasets{cal_dat}.ft_raw_no_offset,1)
        r_ft = datasets{cal_dat}.ft_raw_no_offset(smpl,:);
        pi_known_ft = added_mass(cal_dat);
        g = datasets{cal_dat}.acc(smpl,:);
        regr_ft = [ kron(r_ft,eye(3,3)) -kron(g,eye(3,6))*static_force_regressor ];
        kt_ft = static_force_matrix(pi_known_ft)*g';
        acc_ft = acc_ft + regr_ft'*kt_ft;
        cov_ft = cov_ft + regr_ft'*regr_ft;
    end
end

x_ft = pinv(cov_ft)*acc_ft;
C_ft = reshape(x_ft(1:18),3,6);
unknown_mass = x_ft(19);

%%
% Check predicted force on calibration datasets 
figure
for i=1:n_datasets
    fprintf(['Fitting  force for dataset ' datasets{i}.name '\n']);
    datasets{i}.predicted_force = (C_ft*datasets{i}.ft_raw_no_offset')';
    [datasets{i}.p_force_predicted,datasets{i}.predicted_force_refitted]   = ellipsoidfit_smart(datasets{i}.predicted_force(1:subsampling:end,:),datasets{i}.acc(1:subsampling:end,:));
    subplot(square_side,square_side,i);
    plot_ellipsoid_im(datasets{i}.p_force_predicted);
    plot3_matrix(datasets{i}.predicted_force(1:subsampling:end,:))
    plot_ellipsoid_im(datasets{i}.p_force);
    plot3_matrix(datasets{i}.ft_raw_no_offset(1:subsampling:end,1:3))
    axis equal
    title('Ellipsoid fitted for predicted forces');
end

%% Validation inertial parameters estimation
for i=1:n_datasets
    fprintf(strcat('Estimating added mass for dataset  ',datasets_names{i},'\n'));
    % estimate inertial parameters for this dataset
    cov_par = zeros(1,1);
    kt_par = zeros(1,1);
    for smpl = 1:size(datasets{i}.predicted_force,1)
        mass_regressor = datasets{i}.acc(smpl,:)';
        cov_par = cov_par + mass_regressor'*mass_regressor;
        kt_par = kt_par + mass_regressor'*datasets{i}.predicted_force(smpl,:)';
    end
    cov_var_inv = pinv(cov_par);
    datasets{i}.mass_estimated{j} = cov_var_inv*kt_par;
end

%% Validation inertial parameters estimation
for i=1:n_datasets
    fprintf(strcat('Estimating added mass for dataset  ',datasets_names{i},' with old calibration matrix\n'));
    % estimate inertial parameters for this dataset
    cov_par = zeros(1,1);
    kt_par = zeros(1,1);
    for smpl = 1:size(datasets{i}.predicted_force,1)
        mass_regressor = datasets{i}.acc(smpl,:)';
        cov_par = cov_par + mass_regressor'*mass_regressor;
        kt_par = kt_par + mass_regressor'*datasets{i}.ft_raw_no_offset(smpl,1:3)';
    end
    cov_var_inv = pinv(cov_par);
    datasets{i}.mass_estimated_old{j} = cov_var_inv*kt_par;
end

end

%% Plot mass estimation results 
added_masses_true_values = zeros(1,n_datasets);
added_masses_estimated_via_old_calibration_matrix = zeros(1,n_datasets);
added_masses_estimated_via_new_calibration_matrix = {};
for j=2:n_datasets
    added_masses_estimated_via_new_calibration_matrix{j} = zeros(1,n_datasets);
end

%% this is working just because the first dataset has zero added mass
for i=1:n_datasets
    added_masses_true_values(i) = added_mass(i);
    added_masses_estimated_via_old_calibration_matrix(i) = datasets{i}.mass_estimated_old{2}-datasets{1}.mass_estimated_old{2};
    for j=2:n_datasets
        added_masses_estimated_via_new_calibration_matrix{j}(i) = ...
            datasets{i}.mass_estimated{j}-datasets{1}.mass_estimated{j};
    end
end

figure
scatter(1:n_datasets,added_masses_true_values,50,'green');
hold on;
scatter(1:n_datasets,added_masses_estimated_via_old_calibration_matrix,50,'red');
hold on;
cmap = colormap;
for j=2:n_datasets
    scatter(1:n_datasets,added_masses_estimated_via_new_calibration_matrix{j},50,cmap(5*j,1:3));
    hold on;
end


figure
scatter(1:n_datasets,abs(added_masses_estimated_via_old_calibration_matrix-added_masses_true_values),50,'red');
hold on
cmap = colormap;
for j=2:n_datasets
    scatter(1:n_datasets,abs(added_masses_estimated_via_new_calibration_matrix{j}-added_masses_true_values),50,cmap(5*j,1:3));
end





% 
% %% Validation plots
% 
% for val_dat = validation_datasets
% %     for smpl = 1:size(datasets{cal_dat}.foot_no_offset,1)
% %     end
%     datasets{val_dat}.calibrated_leg_no_offset = datasets{val_dat}.leg_no_offset*C_leg';
%     datasets{val_dat}.calibrated_leg_no_offset_original  = datasets{val_dat}.leg_no_offset*original_calibration_matrix_leg';
% end
% 
% light_red   = [255.0 192.0 203.0]./255.0;
% dark_red = [255.0 0 0]./255.0;
% light_green = [144 238.0 144]./255.0;
% dark_green  = [0 128.0 0]./255.0;
% 
% figure
% for val_dat = 5:n_datasets
%     fprintf(strcat('Plotting calibrated ellipsoids for foot force, dataset',datasets_names{val_dat}));
% 
% %   for smpl = 1:size(datasets{cal_dat}.foot_no_offset,1)
% %   end
%     hold on;
%     subplot(2,2,val_dat-4);
%     p_calib = ellipsoidfit_smart(datasets{val_dat}.calibrated_foot_no_offset(1:subsampling:end,1:3),datasets{val_dat}.acc(1:subsampling:end,:));
%     ellipsoid_ecc(p_calib)
%     hold on
%     hold on
%     hold on
%     p_calib_original = ellipsoidfit_smart(datasets{val_dat}.calibrated_foot_no_offset_original(1:subsampling:end,1:3),datasets{val_dat}.acc(1:subsampling:end,:));
%     ellipsoid_ecc(p_calib_original)
%     plot_ellipsoid_im_color(p_calib_original,light_red);
%     plot_ellipsoid_im_color(p_calib,light_green);
%     plot3_matrix_color(datasets{val_dat}.calibrated_foot_no_offset_original(:,1:3),dark_red);
%     plot3_matrix_color(datasets{val_dat}.calibrated_foot_no_offset(:,1:3),dark_green);
% 
%     view([45 10]);
%     
%     axis equal;
%     axis([-11 11 -11 11 -11 11]);
%     xlabel('f_x (N)')
%     ylabel('f_y (N)')
%     zlabel('f_z (N)')
% end
% 
% figure
% for val_dat = 5:n_datasets
%     fprintf(strcat('Plotting calibrated ellipsoids for leg force',datasets_names{val_dat}));
% %     for smpl = 1:size(datasets{cal_dat}.foot_no_offset,1)
% %     end
%     hold on;
%     subplot(2,2,val_dat-4);
%     p_calib_original = ellipsoidfit_smart(datasets{val_dat}.calibrated_leg_no_offset_original(1:subsampling:end,1:3),datasets{val_dat}.acc(1:subsampling:end,:));
%       ellipsoid_ecc(p_calib_original)
%     p_calib = ellipsoidfit_smart(datasets{val_dat}.calibrated_leg_no_offset(1:subsampling:end,1:3),datasets{val_dat}.acc(1:subsampling:end,:));
% 
%     ellipsoid_ecc(p_calib)
%     plot_ellipsoid_im_color(p_calib_original,light_red);
%     plot3_matrix_color(datasets{val_dat}.calibrated_leg_no_offset_original(:,1:3),dark_red);
% %    plot_ellipsoid_im_color(p_calib,light_green);
% %    plot3_matrix_color(datasets{val_dat}.calibrated_leg_no_offset(:,1:3),dark_green);
%     view([45 10]);
%     axis([-40 40 -60 60 -40 40]);
%     axis equal
%     xlabel('f_x (N)')
%     ylabel('f_y (N)')
%     zlabel('f_z (N)')
% end
% 
% %% Validation inertial parameters estimation
% 
% 
% for i=1:n_datasets
%     fprintf(strcat('Estimating inertial parameters for dataset ',datasets_names{i}),'\n');
%     % estimate inertial parameters for this dataset
%     cov_par = zeros(4,4);
%     kt_par_leg = zeros(4,1);
%     kt_par_leg_old = zeros(4,1);
%     for smpl = 1:size(datasets{i}.calibrated_leg_no_offset,1)
%         inertial_param_regressor = inertial_parameter_static_regressor(datasets{i}.acc(smpl,:)');
%         cov_par = cov_par + inertial_param_regressor'*inertial_param_regressor;
%         kt_par_leg = kt_par_leg + inertial_param_regressor'*datasets{i}.calibrated_leg_no_offset(smpl,:)';
%         kt_par_leg_old = kt_par_leg_old + inertial_param_regressor'*datasets{i}.calibrated_leg_no_offset_original(smpl,:)';
%     end
%     cov_var_inv = pinv(cov_par);
%     datasets{i}.pi_leg = cov_var_inv*kt_par_leg;
%     datasets{i}.pi_leg_old = cov_var_inv*kt_par_leg_old;    
% end
% 
% for i=1:n_datasets
%     datasets{i}.pi_leg_added_known          = [added_mass(i) added_first_moment_of_mass_wrt_leg{i}];
%     datasets{i}.pi_leg_added_estimated      = datasets{i}.pi_leg      - pi_unknown_leg;
%     datasets{i}.pi_leg_added_estimated_old  = datasets{i}.pi_leg_old  - pi_unknown_leg;
% end
% 
% figure
% for i=1:n_datasets
%     subplot(2,4,i);
%     x = [1:3,6:8,11:13,16:18];
%     y = [datasets{i}.pi_leg_added_known(1),datasets{i}.pi_leg_added_estimated(1),datasets{i}.pi_leg_added_estimated_old(1),...
%          datasets{i}.pi_leg_added_known(2),datasets{i}.pi_leg_added_estimated(2),datasets{i}.pi_leg_added_estimated_old(2),...
%          datasets{i}.pi_leg_added_known(3),datasets{i}.pi_leg_added_estimated(3),datasets{i}.pi_leg_added_estimated_old(3),...
%          datasets{i}.pi_leg_added_known(4),datasets{i}.pi_leg_added_estimated(4),datasets{i}.pi_leg_added_estimated_old(4)];
%      bar(x,y);
% end
% 
% % %% Diagnostics plots for debugging the local drift of the offset
% % figure
% % for i = 1:n_datasets
% %     subplot(square_side,square_side,i);
% %     plot3_matrix(datasets{i}.ft_leg_projected_norm(1:subsampling:end,:));
% %     hold on
% %     scatter3_matrix(datasets{i}.ft_leg_projected_norm(datasets{i}.return_point==1.0,:));
% %     axis equal
% % end
% % 
% % figure
% % for i = 1:n_datasets
% %     subplot(square_side,square_side,i);
% %     plot3_matrix(datasets{i}.ft_foot_projected_norm(1:subsampling:end,:));
% %     hold on
% %     scatter3_matrix(datasets{i}.ft_foot_projected_norm(datasets{i}.return_point==1.0,:));
% %     axis equal
% % end
% % 
% % figure
% % for i = 1:n_datasets
% %     subplot(square_side,square_side,i);
% %     plot3_matrix(datasets{i}.acc(1:subsampling:end,:));
% %     hold on
% %     scatter3_matrix(datasets{i}.acc(datasets{i}.return_point==1.0,:));
% %     axis equal
% % end

 dlmwrite('validation_force_matrix.csv',C_ft,',')
