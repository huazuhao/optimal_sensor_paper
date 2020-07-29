clc;
clear;

rng(4);



% Problem parameters
n     = 100;  % Number of interior spatial dofs
N     = 100;  % Number of temporal dofs
gamma = 0.5; % Robin coefficient, this is part of the boundary condition

[M, A, B, R] = buildFEM(n, gamma);


total_a = full(-1*A);
total_m = full(M);


%define a new basis that we want to project to
new_basis_rank = 6;
selected_j = [0,1,2,1,2,2];
selected_k = [0,1,2,0,0,1];
new_basis_vectors = [];

for index = 1:new_basis_rank
    
    j_index = selected_j(index);
    k_index = selected_k(index);
    
    a_vector = haar_hat_function( j_index,k_index,n );
    
    new_basis_vectors = [new_basis_vectors,a_vector];
    
end

new_basis = inv(total_m)*new_basis_vectors;



%we now begin to select sensors
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
            projected_og = new_basis'*og_matlab*new_basis;
            [u_2,s_2,v_2] = svd(projected_og);
            objective_value = sum(log(diag(s_2)));
            log_det_og_matlab(index) = objective_value;
            
        end
    end
    
    [~,optimal_sensor_placement] = max(log_det_og_matlab);
    selected_obs = [selected_obs,optimal_sensor_placement]
    
end

projection_matrix = new_basis*inv(new_basis'*new_basis)*new_basis';
test_projection_matrix = projection_matrix-projection_matrix*projection_matrix;




all_index = zeros(1,n+2);
all_index(selected_obs) = 1;
figure;

scatter(selected_obs,all_index(selected_obs),50,'d','filled');


hold on;
plot(new_basis)
xlabel('position')
ylabel('haar function value')
legend('sensor locations','haar basis');

save('exact_method_6_choose_5.mat','all_index');
haar_index = [selected_j;selected_k];
save('exact_haar_index_6.mat','haar_index');