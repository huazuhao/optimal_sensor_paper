function [ integral_result ] = haar_hat_function( haar_j,haar_k,num_of_interior_points )
%HAAR_HAT_FUNCTION Summary of this function goes here
%   Detailed explanation goes here

    haar_one_begin = haar_k/(2^haar_j);
    haar_one_end = (haar_k+1/2)/(2^haar_j);
    haar_minus_one_begin = (haar_k+1/2)/(2^haar_j);
    haar_minus_one_end = (haar_k+1)/(2^haar_j);
    haar_factor = 2^(haar_j/2);

    
    
    num_of_divided_interval = num_of_interior_points+1; 

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

    %now, do the intergral
    integral_result = [];
    for index = 1:num_of_interior_points

        integral_result_case_one = 0;
        integral_result_case_two = 0;
        integral_result_case_three = 0;
        integral_result_case_four = 0;


        left = left_end_hat(index);
        middle = hat_middle_point(index);
        right = right_end_hat(index);

        pos_left_hat_int_begin = max(haar_one_begin,left);
        pos_left_hat_int_end = min(haar_one_end,middle);
        if pos_left_hat_int_begin>=left && pos_left_hat_int_begin<=middle
            if pos_left_hat_int_end>=left && pos_left_hat_int_end<=middle
                int_func = @(x) haar_factor*(x-left)/(1/num_of_divided_interval);
                integral_result_case_one = integral(int_func,pos_left_hat_int_begin,pos_left_hat_int_end);
            end
        end

        neg_left_hat_int_begin = max(haar_minus_one_begin,left);
        neg_left_hat_int_end = min(haar_minus_one_end,middle);
        if neg_left_hat_int_begin>=left && neg_left_hat_int_begin<=middle
            if neg_left_hat_int_end>=left && neg_left_hat_int_end<=middle
                int_func = @(x) haar_factor*-1*(x-left)/(1/num_of_divided_interval);
                integral_result_case_two = integral(int_func,neg_left_hat_int_begin,neg_left_hat_int_end);
            end
        end

        pos_right_hat_int_begin = max(haar_one_begin,middle);
        pos_right_hat_int_end = min(haar_one_end,right);
        if pos_right_hat_int_begin>=middle && pos_right_hat_int_begin<=right
            if pos_right_hat_int_end>=middle && pos_right_hat_int_end<=right
                int_func = @(x) haar_factor*(right-x)/(1/num_of_divided_interval);
                integral_result_case_three = integral(int_func,pos_right_hat_int_begin,pos_right_hat_int_end);
            end
        end

        neg_right_hat_int_begin = max(haar_minus_one_begin,middle);
        neg_right_hat_int_end = min(haar_minus_one_end,right);
        if neg_right_hat_int_begin>=middle && neg_right_hat_int_begin<=right
            if neg_right_hat_int_end>=middle && neg_right_hat_int_end<=right
                int_func = @(x) haar_factor*-1*(right-x)/(1/num_of_divided_interval);
                integral_result_case_four = integral(int_func,neg_right_hat_int_begin,neg_right_hat_int_end);
            end
        end

        total_integral = integral_result_case_one+integral_result_case_two+...
                         integral_result_case_three+integral_result_case_four;

        integral_result = [integral_result;total_integral];

    end
    
    
    
    %left end
    integral_result_one = 0;
    integral_result_two = 0;
    pos_hat_int_begin = max(haar_one_begin,0);
    pos_hat_int_end = min(haar_one_end,1/num_of_divided_interval);
    if pos_hat_int_begin>=0 && pos_hat_int_begin<=1/num_of_divided_interval
        if pos_hat_int_end>=0 && pos_hat_int_end<=1/num_of_divided_interval
            int_func = @(x) haar_factor*(1/num_of_divided_interval-x)/(1/num_of_divided_interval);
            integral_result_one = integral(int_func,pos_hat_int_begin,pos_hat_int_end);
        end
    end
    neg_hat_int_begin = max(haar_minus_one_begin,0);
    neg_hat_int_end = min(haar_minus_one_end,1/num_of_divided_interval);
    if neg_hat_int_begin>=0 && neg_hat_int_begin<=1/num_of_divided_interval
        if neg_hat_int_end>=0 && neg_hat_int_end<=1/num_of_divided_interval
            int_func = @(x) -1*haar_factor*(1/num_of_divided_interval-x)/(1/num_of_divided_interval);
            integral_result_two = integral(int_func,neg_hat_int_begin,neg_hat_int_end);
        end
    end
    total_integral = integral_result_one+integral_result_two;
    integral_result = [total_integral;integral_result];

    %right end
    integral_result_one = 0;
    integral_result_two = 0;
    pos_hat_int_begin = max(haar_one_begin,1-1/num_of_divided_interval);
    pos_hat_int_end = min(haar_one_end,1);
    if pos_hat_int_begin>=1-1/num_of_divided_interval && pos_hat_int_begin<=1
        if pos_hat_int_end>=1-1/num_of_divided_interval && pos_hat_int_end<=1
            int_func = @(x) haar_factor*(x-(1-1/num_of_divided_interval))/(1/num_of_divided_interval);
            integral_result_one = integral(int_func,pos_hat_int_begin,pos_hat_int_end);
        end
    end
    neg_hat_int_begin = max(haar_minus_one_begin,1-1/num_of_divided_interval);
    neg_hat_int_end = min(haar_minus_one_end,1);
    if neg_hat_int_begin>=1-1/num_of_divided_interval && neg_hat_int_begin<=1
        if neg_hat_int_end>=1-1/num_of_divided_interval && neg_hat_int_end<=1
            int_func = @(x) -1*haar_factor*(x-(1-1/num_of_divided_interval))/(1/num_of_divided_interval);
            integral_result_two = integral(int_func,neg_hat_int_begin,neg_hat_int_end);
        end
    end
    total_integral = integral_result_one+integral_result_two;
    integral_result = [integral_result;total_integral];
    


end

