clc;
clear;

rng(2);



% Problem parameters
n     = 100;  % Number of interior spatial dofs
N     = 100;  % Number of temporal dofs
gamma = 0.5; % Robin coefficient, this is part of the boundary condition
rank_cut_off_for_integration = 30;


[M, A, B, R] = buildFEM(n, gamma);


total_a = full(-1*A);
total_m = full(M);


[w,l] = get_w_l_matrix_drew(M,A);
wt = w';


%define a new basis that we want to project to
random_matrix = rand(n+2,n+2);
[new_basis,r] = qr(random_matrix);
new_basis_rank = 5;
new_basis = new_basis(:,1:new_basis_rank);


special_w = new_basis'*w;
special_w = special_w(:,n+3-rank_cut_off_for_integration:n+2);



selected_obs = [];
%now, I initialize the u,s,v matrix by doing a svd on a zero matrix
[u_old,s_old,v_old] = svd(zeros(new_basis_rank,new_basis_rank));

if length(selected_obs)>0
    c_matrix = zeros(n+2,n+2);
    for index = 1:length(selected_obs)
        already_selected = selected_obs(index);
        c_matrix(already_selected,already_selected) = 1;
    end
    %compute the svd of v^tOGv
    og_matlab = lyap(transpose(total_a),c_matrix'*c_matrix,[],transpose(total_m));
    projected_og = new_basis'*og_matlab*new_basis;
    [u_old,s_old,v_old] = svd(projected_og);
end


%%
for greedy_iteration = 1:10
    
    disp(greedy_iteration)
    
    %check where we should pick next
    next_potential_placement = ones(1,n+2);
    for index = 1:length(selected_obs)
        already_selected = selected_obs(index);
        next_potential_placement(already_selected) = 0;
    end
    
    log_det_og = ones(1,n+2)*-1e9;
    log_det_og_matlab = ones(1,n+2)*-1e9;
    
    p = (eye(new_basis_rank)-u_old*u_old')*special_w;
    [p,r] = qr(p,0);
    up = [u_old,p];
    r_size = size(r);
    ua_size = size(u_old'*special_w);
    k_left = [eye(new_basis_rank),u_old'*special_w;zeros(r_size(1),ua_size(1)),r];

    q = (eye(new_basis_rank)-v_old*v_old')*special_w;
    [q,r] = qr(q,0);
    vq = [v_old,q];
    r_size = size(r);
    vb_size = size(v_old'*special_w);
    k_right = [eye(new_basis_rank),v_old'*special_w;zeros(r_size(1),vb_size(1)),r];

    s_old_size = size(s_old);
    
    c_old = zeros(n+2,n+2);
    for index = 1:length(selected_obs)
        already_selected = selected_obs(index);
        c_old(already_selected,already_selected) = 1;
    end
    
    for index = 1:n+2
        
        if next_potential_placement(index) ~= 0
        
             %the first step is to get the wtctcw_diff matrix correct
%             c_new = zeros(n+2,n+2);
%             c_new(index,index) = 1;
%             newnew = c_new'*total_m*c_new;
%             newold = c_new'*total_m*c_old;
%             oldnew = c_old'*total_m*c_new;
%             middle_diff = wt*(newnew+newold+oldnew)*w;
%             
            [newnew,newold,oldnew] = compute_middle_diff(index,selected_obs,wt,w,total_m,n,rank_cut_off_for_integration);
            middle_diff = newnew+newold+oldnew;
            
            %the second step is to do the integration
            integration_diff = zeros(n+2,n+2);
            for row = n+3-rank_cut_off_for_integration:n+2
                for column = n+3-rank_cut_off_for_integration:n+2 
                    each_entry = middle_diff(row,column)*-1/(l(row,row)+l(column,column));
                    integration_diff(row,column) = each_entry;
                end
            end
            in_term = integration_diff(n+3-rank_cut_off_for_integration:n+2,...
                                n+3-rank_cut_off_for_integration:n+2);

            in_term_size = size(in_term);


            %form the middle k matrix on matthew brand's paper
            k = k_left*[s_old,zeros(s_old_size(1),in_term_size(1));zeros(in_term_size(1),s_old_size(1)),in_term]*k_right';


            %we perform the svd of the k matrix for computing the objective
            [u_2,s_2,v_2] = svd(k);

            s_temp = s_2(1:new_basis_rank,1:new_basis_rank);
            objective_value = sum(log(diag(s_temp)));
            %objective_value = sum(log(diag(s_2)));
            log_det_og(index) = objective_value;
            
            
            %let's implement a matlab call to check the effect of adding
            %this new sensor
%             c_matrix = zeros(n+2,n+2);
%             for index2 = 1:length(selected_obs)
%                 already_selected = selected_obs(index2);
%                 c_matrix(already_selected,already_selected) = 1;
%             end
%             c_matrix(index,index) = 1;
%             og_matlab = lyap(transpose(total_a),c_matrix'*total_m*c_matrix,[],transpose(total_m));
%             projected_og = new_basis'*og_matlab*new_basis;
%             [u_2,s_2,v_2] = svd(projected_og);
%             objective_value = sum(log(diag(s_2)));
%             log_det_og_matlab(index) = objective_value;
            

        end%end of the if statement
    end% end of the for loop for placing new sensors
     
%     select_logical = next_potential_placement>0;
%     cleaned_log_det_og = log_det_og(:,select_logical);
%     cleaned_log_det_og_matlab = log_det_og_matlab(:,select_logical);
%     %plot(cleaned_log_det_og)
%     %hold on;
%     %plot(cleaned_log_det_og_matlab)
%     cleaned_log_det_og-cleaned_log_det_og_matlab
    
    [~,optimal_sensor_placement] = max(log_det_og);
    selected_obs = [selected_obs,optimal_sensor_placement]
    
    c_new = zeros(n+2,n+2);
    c_new(optimal_sensor_placement,optimal_sensor_placement) = 1;
    newnew = c_new'*total_m*c_new;
    newold = c_new'*total_m*c_old;
    oldnew = c_old'*total_m*c_new;
    
    middle_diff = w'*(newnew+newold+oldnew)*w;
    
    for row = n+3-rank_cut_off_for_integration:n+2
        for column = n+3-rank_cut_off_for_integration:n+2
            entry = middle_diff(row,column)*-1/(l(row,row)+l(column,column));
            in_term(row,column) = entry;
        end
    end
    in_term = in_term(n+3-rank_cut_off_for_integration:n+2,n+3-rank_cut_off_for_integration:n+2);
    in_term_size = size(in_term);
    k = k_left*[s_old,zeros(s_old_size(1),in_term_size(1));zeros(in_term_size(1),s_old_size(1)),in_term]*k_right';
    
    %now, we do the svd of the k matrix
    [u_2,s_2,v_2] = svd(k);
    %compute_the_objective
    u_update_once = up*u_2;
    s_update_once = s_2;
    v_update_once = vq*v_2;
    
    u_old = u_update_once(:,1:new_basis_rank);
    s_old = s_update_once(1:new_basis_rank,1:new_basis_rank);
    v_old = v_update_once(:,1:new_basis_rank);
    
end



all_index = zeros(1,n+2);
all_index(selected_obs) = 1;
figure;
plot(all_index);


save('approximation_method.mat','all_index');