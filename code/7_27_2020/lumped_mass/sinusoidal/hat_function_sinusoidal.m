function [integral_result] = hat_function_sinusoidal(wavenumber_input,num_of_interior_points_input)

    wavenumber = wavenumber_input;
    num_of_interior_points = num_of_interior_points_input;

    num_of_divided_interval = num_of_interior_points+1; 


    %generate left,middle,right points for hat functions
    hat_middle_point = [];
    for index = 1:num_of_interior_points
        point = index*1/num_of_divided_interval;
        hat_middle_point = [hat_middle_point,point];
    end

    left_end_hat = [];
    right_end_hat = [];
    for index = 1:num_of_interior_points
        middle = hat_middle_point(index);
        left = middle-1/num_of_divided_interval;
        right = middle+1/num_of_divided_interval;
        left_end_hat = [left_end_hat,left];
        right_end_hat = [right_end_hat,right];
    end

    %now we do the integral
    integral_result = [];
    for index = 1:num_of_interior_points

        left = left_end_hat(index);
        middle = hat_middle_point(index);
        right = right_end_hat(index);

        %left hat
        int_func = @(x) sqrt(2)*sin(2*pi*wavenumber*x).*(x-left)/(1/num_of_divided_interval);
        left_hat_int = integral(int_func,left,middle);

        %right hat
        int_func = @(x) sqrt(2)*sin(2*pi*wavenumber*x).*(right-x)/(1/num_of_divided_interval);
        right_hat_int = integral(int_func,middle,right);

        total_integral = left_hat_int+right_hat_int;

        integral_result = [integral_result;total_integral];

    end
    
    %left hat basis
    int_func = @(x) sqrt(2)*sin(2*pi*wavenumber*x).*(1/num_of_divided_interval-x)/(1/num_of_divided_interval);
    left_hat_int = integral(int_func,0,1/num_of_divided_interval);
    
    integral_result = [left_hat_int;integral_result];
    
    %right hat basis
    int_func = @(x) sqrt(2)*sin(2*pi*wavenumber*x).*(x-(1-1/num_of_divided_interval))/(1/num_of_divided_interval);
    right_hat_int = integral(int_func,1-1/num_of_divided_interval,1);
    
    integral_result = [integral_result;right_hat_int];
    
    
    
end
    