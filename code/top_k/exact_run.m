clc;
clear;

rng(0);



% Problem parameters
n     = 100;  % Number of interior spatial dofs
N     = 100;  % Number of temporal dofs
gamma = 0.5; % Robin coefficient, this is part of the boundary condition
rank_cut_off = 10;

[M, A, B, R] = buildFEM(n, gamma);

total_a = full(-1*A);
total_m = full(M);

selected_obs = [];

%%
for greedy_iteration = 1:6

    disp(greedy_iteration)
    
    %check where we should pick next
    next_potential_placement = ones(1,n+2);
    for index = 1:length(selected_obs)
        already_selected = selected_obs(index);
        next_potential_placement(already_selected) = 0;
    end
    
    log_det_og_matlab = ones(1,n+2)*(-1*1e9);
    %log_det_og_matlab = zeros(1,n+2);
    
    for index = 1:n+2
    
        if next_potential_placement(index) ~= 0
            
            new_c_matrix = zeros(n+2,n+2);
            for index_2 = 1:length(selected_obs)
                non_zero_index = selected_obs(index_2);
                new_c_matrix(non_zero_index, non_zero_index) = 1;
            end
            new_c_matrix(index,index) = 1;
            
            og_matlab = lyap(transpose(total_a),new_c_matrix'*total_m*new_c_matrix,[],transpose(total_m));
            %og_matlab = lyap(transpose(total_a),new_c_matrix'*new_c_matrix,[],transpose(total_m));
            
            [~,s_matlab,~] = svd(og_matlab);
            s_low_rank = s_matlab(1:rank_cut_off,1:rank_cut_off);
            log_det_low_rank = sum(log(diag(s_low_rank)));
            log_det_og_matlab(index) = log_det_low_rank;
            
        end %finish the if statement
    end%finish checking every possible sensor placement for one loop of greedy algorithm
    
    [~,optimal_sensor_placement] = max(log_det_og_matlab);
    selected_obs = [selected_obs,optimal_sensor_placement];
    
end

all_index = zeros(1,n+2);
all_index(selected_obs) = 1;
        
figure
plot(all_index);
xlim([1 length(all_index)]);

save('exact_method.mat','all_index');

