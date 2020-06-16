clc;
clear;

rng(0);


% Problem parameters
n     = 300;  % Number of interior spatial dofs
N     = 300;  % Number of temporal dofs
gamma = 0.5; % Robin coefficient, this is part of the boundary condition
rank_cut_off_for_integration = 100;
rank_cut_off = 10;

[M, A, B, R] = buildFEM(n, gamma);

total_a = full(-1*A);
total_m = full(M);

[w,l] = get_w_l_matrix_drew(M,A);
wt = w';
special_w = w(:,1:rank_cut_off_for_integration);


selected_obs = [];
%now, I initialize the u,s,v matrix by doing a svd on a zero matrix
[u_old,s_old,v_old] = svd(zeros(n+2,n+2));

u_old = u_old(:,1:rank_cut_off_for_integration);
s_old = s_old(1:rank_cut_off_for_integration,1:rank_cut_off_for_integration);
v_old = v_old(:,1:rank_cut_off_for_integration);



for greedy_iteration = 1:10
    
    disp(greedy_iteration)
    
    
    %check where we should pick next
    next_potential_placement = ones(1,n+2);
    for index = 1:length(selected_obs)
        already_selected = selected_obs(index);
        next_potential_placement(already_selected) = 0;
    end
    
    %some preparation for the next sensor
    log_det_og = ones(1,n+2)*-1e9;
    log_det_og_matlab = ones(1,n+2)*-1e9;
    
    p = (eye(n+2)-u_old*u_old')*special_w;
    [p,r] = qr(p,0);
    up = [u_old,p];
    r_size = size(r);
    ua_size = size(u_old'*special_w);
    k_left = [eye(rank_cut_off_for_integration),u_old'*special_w;zeros(r_size(1),ua_size(1)),r];

    q = (eye(n+2)-v_old*v_old')*special_w;
    [q,r] = qr(q,0);
    vq = [v_old,q];
    r_size = size(r);
    vb_size = size(v_old'*special_w);
    k_right = [eye(rank_cut_off_for_integration),v_old'*special_w;zeros(r_size(1),vb_size(1)),r];
    
    s_old_size = size(s_old);
    
    c_old = zeros(n+2,n+2);
    for index = 1:length(selected_obs)
        already_selected = selected_obs(index);
        c_old(already_selected,already_selected) = 1;
    end
    
        for index = 1:n+2
        
            if next_potential_placement(index) ~= 0

                c_new = zeros(n+2,n+2);
                c_new(index,index) = 1;
                newnew = c_new'*total_m*c_new;
                newold = c_new'*total_m*c_old;
                oldnew = c_old'*total_m*c_new;

                middle_diff = w'*(newnew+newold+oldnew)*w;

                in_term = zeros(n+2,n+2);
                for row = 1:rank_cut_off_for_integration
                    for column = 1:rank_cut_off_for_integration
                        entry = middle_diff(row,column)*-1/(l(row,row)+l(column,column));
                        in_term(row,column) = entry;
                    end
                end
                in_term = in_term(1:rank_cut_off_for_integration,1:rank_cut_off_for_integration);


                s_old_size = size(s_old);
                in_term_size = size(in_term);

                k = k_left*[s_old,zeros(s_old_size(1),in_term_size(1));zeros(in_term_size(1),s_old_size(1)),in_term]*k_right';


                %now, we do the svd of the k matrix
                [u_2,s_2,v_2] = svd(k);

                s_temp = s_2(1:rank_cut_off,1:rank_cut_off);
                objective_value = sum(log(diag(s_temp)));
                log_det_og(index) = objective_value;
                
                %now, let's compute the objective with the matlab call
%                 c_ab_initio = zeros(n+2,n+2);
%                 for index2 = 1:length(selected_obs)
%                     already_selected = selected_obs(index2);
%                     c_ab_initio(already_selected,already_selected) = 1;
%                 end
%                 c_ab_initio(index,index) = 1;
%                 og_matlab = lyap(transpose(total_a),c_ab_initio'*total_m*c_ab_initio,[],transpose(total_m));
% 
%                 og_matlab-up*k*vq';
% 
%                 [u_2_matlab,s_2_matlab,v_2_matlab] = svd(og_matlab);
%                 s_temp = s_2_matlab(1:rank_cut_off,1:rank_cut_off);
%                 objective_value = sum(log(diag(s_temp)));
%                 log_det_og_matlab(index) = objective_value;


             end %end if statement
             
        end %end new sensor for loop
        
%     select_logical = next_potential_placement>0;
%     cleaned_log_det_og = log_det_og(:,select_logical);
%     cleaned_log_det_og_matlab = log_det_og_matlab(:,select_logical);
%     cleaned_log_det_og-cleaned_log_det_og_matlab
    
        
    [~,optimal_sensor_placement] = max(log_det_og);
    selected_obs = [selected_obs,optimal_sensor_placement];
    
    c_new = zeros(n+2,n+2);
    c_new(optimal_sensor_placement,optimal_sensor_placement) = 1;
    newnew = c_new'*total_m*c_new;
    newold = c_new'*total_m*c_old;
    oldnew = c_old'*total_m*c_new;

    middle_diff = w'*(newnew+newold+oldnew)*w;
    
    for row = 1:rank_cut_off_for_integration
        for column = 1:rank_cut_off_for_integration
            entry = middle_diff(row,column)*-1/(l(row,row)+l(column,column));
            in_term(row,column) = entry;
        end
    end
    in_term = in_term(1:rank_cut_off_for_integration,1:rank_cut_off_for_integration);
    in_term_size = size(in_term);
    k = k_left*[s_old,zeros(s_old_size(1),in_term_size(1));zeros(in_term_size(1),s_old_size(1)),in_term]*k_right';
            
    %now, we do the svd of the k matrix
    [u_2,s_2,v_2] = svd(k);
    %compute_the_objective
    u_update_once = up*u_2;
    s_update_once = s_2;
    v_update_once = vq*v_2;
    
    u_old = u_update_once(:,1:rank_cut_off_for_integration);
    s_old = s_update_once(1:rank_cut_off_for_integration,1:rank_cut_off_for_integration);
    v_old = v_update_once(:,1:rank_cut_off_for_integration);
    
    [q_new_u,r_new_u] = qr(u_update_once);
    [q_new_v,r_new_v] = qr(v_update_once);
    s_new = r_new_u*s_2*r_new_v';
    u_old = q_new_u(:,1:rank_cut_off_for_integration);
    v_old = q_new_v(:,1:rank_cut_off_for_integration);
    s_old = s_new(1:rank_cut_off_for_integration,1:rank_cut_off_for_integration);
    
end


all_index = zeros(1,n+2);
all_index(selected_obs) = 1;
        
figure
plot(all_index);
xlim([1 length(all_index)]);

save('approximation_method.mat','all_index');