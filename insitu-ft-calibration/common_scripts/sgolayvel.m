function [ ret_smooth, ret_vel ] = sgolayvel( input_mat, order, window )
%SGOLAYVEL Return a matrix where the columns are the velocity with sgolay
%   of the columns
    [b,g] = sgolay(order,window);
    HalfWin = ((window+1)/2) -1;
    
    n_samples = size(input_mat,1);
    n_channels = size(input_mat,2);
    
    ret_smooth = input_mat;
    ret_vel = input_mat;
       
    for n_ch = 1:n_channels,
        for n_sm = (window+1)/2:(n_samples-(window/2)),
            ret_smooth(n_sm,n_ch) = dot(g(:,1),input_mat(n_sm - HalfWin: n_sm + HalfWin,n_ch));
        end
    end  
    
    for n_ch = 1:n_channels,
        for n_sm = (window+1)/2:(n_samples-(window/2)),
             ret_vel(n_sm,n_ch) = dot(g(:,2),input_mat(n_sm - HalfWin: n_sm + HalfWin,n_ch));
        end
    end
    
    %return only the filtered values
    ret_vel = ret_vel((window+1)/2:(n_samples-(window/2)),:);
    ret_smooth = ret_smooth((window+1)/2:(n_samples-(window/2)),:);

end

