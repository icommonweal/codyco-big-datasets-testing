function [ p_implicit, output_data_refitted ] = ellipsoidfit_smart( output_elliptic_data, input_circular_data )
%ellipsoidfit_smart 
    assert( size(output_elliptic_data,1) == size(input_circular_data,1) );
    n_samples = size(input_circular_data,1);
    assert( size(output_elliptic_data,2) == size(input_circular_data,2) );
    assert( size(output_elliptic_data,2) == 3 );
%   This is how it works in theory, but it takes to much memory
%     known_terms = output_elliptic_data';
%     known_terms = known_terms(:);
%     regressor = [kron(input_circular_data,eye(3)),repmat(eye(3),n_samples)];
    cov = zeros(12,12);
    acc = zeros(12,1);
    
    for i=1:n_samples
        input = input_circular_data(i,:);
        output = output_elliptic_data(i,:);
        regr = [kron(input,eye(3)),eye(3)];
        cov = cov + transpose(regr)*regr;
        acc = acc + transpose(regr)*(output');
    end
    
    x = pinv(cov)*acc;
        
    a = mean(sqrt(sum(input_circular_data.^2,2)));
    C = a*reshape(x(1:9),3,3);
    o = x(10:12);
    
%     figure
%     plot3(input_circular_data(:,1),input_circular_data(:,2),input_circular_data(:,3),'green');
%     hold on
%     plot3(output_elliptic_data(:,1),output_elliptic_data(:,2),output_elliptic_data(:,3),'red');
% 
     output_data_refitted = (C./a*input_circular_data'+o*ones(1,n_samples))';
%     hold on
%     plot3(output_from_input(:,1),output_from_input(:,2),output_from_input(:,3),'blue');
   
    %check ellipsoid_im2ex for the scheme used for implicit expression of
    %ellipsoid

    S = inv(C);
    Q = S'*S;
    v = -o'*Q;
    p_implicit = zeros(10,1);
    p_implicit(1) = Q(1,1);
    p_implicit(2) = Q(2,2);
    p_implicit(3) = Q(3,3);
    p_implicit(4) = 2*Q(1,2);
    p_implicit(5) = 2*Q(1,3);
    p_implicit(6) = 2*Q(2,3);
    p_implicit(7) = 2*v(1);
    p_implicit(8) = 2*v(2);
    p_implicit(9) = 2*v(3);
    p_implicit(10) = o'*Q*o-1;
   
    
%     figure
% %     plot3_matrix(input_circular_data./a);
%     plot3(output_elliptic_data(:,1),output_elliptic_data(:,2),output_elliptic_data(:,3),'red');
%     plot_ellipsoid_im(p_implicit);

end

