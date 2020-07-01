clc;
clear;

rng(2);



% Problem parameters
n     = 120;  % Number of interior spatial dofs
N     = 120;  % Number of temporal dofs
gamma = 0.5; % Robin coefficient, this is part of the boundary condition

[M, A, B, R] = buildFEM(n, gamma);


total_a = full(-1*A);
total_m = full(M);

selected_obs = [];


for greedy_iteration = 1:5
    
    disp(greedy_iteration)
    
    %check where we should pick next
    next_potential_placement = ones(1,n+2);
    for index = 1:length(selected_obs)
        already_selected = selected_obs(index);
        next_potential_placement(already_selected) = 0;
    end

    log_det_og_matlab = ones(1,n+2)*-1e9;
   
    for index = 1:n+2
        
        if next_potential_placement(index) ~= 0
    
            c_matrix = zeros(n+2,n+2);
            for index2 = 1:length(selected_obs)
                already_selected = selected_obs(index2);
                c_matrix(already_selected,already_selected) = 1;
            end
            c_matrix(index,index) = 1; %effect of new sensor
            
            og_matlab = lyap(transpose(total_a),c_matrix'*total_m*c_matrix,[],transpose(total_m));
            [u_2,s_2,v_2] = svd(og_matlab);
            objective_value = sum(log(diag(s_2)));
            log_det_og_matlab(index) = objective_value;
            
        end
    end
    
    [~,optimal_sensor_placement] = max(log_det_og_matlab);
    selected_obs = [selected_obs,optimal_sensor_placement]
    
end

all_index = zeros(1,n+2);
all_index(selected_obs) = 1;
figure;
plot(all_index);
    
    