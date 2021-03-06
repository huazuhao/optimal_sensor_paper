function [outputArg1,outputArg2] = infinity_series_diff_function(space_dimension_size,mu,max_term)




% Problem parameters
n     = space_dimension_size;  % Number of interior spatial dofs
N     = 10;  % Number of temporal dofs
theta = 0.5;  % Trapezoidal rule, this parameter is used for discretization
gamma = 0.5; % Robin coefficient, this is part of the boundry condition
alpha = 0.0001; % Control penalty parameter, this parameter originates from the objective function
heat_constant = 1; %The constant in the heat equation. 

[M, A, B, R] = buildFEM(n, gamma);
A = heat_constant*A;
B = heat_constant*B;
total_a = full(-1*A);
total_m = full(M);


%given the old sensor and compute the observability gramian
old_sensors = [1,5,9];
old_c_matrix = zeros(n+2,n+2);
for index = 1:length(old_sensors)
    sensor_location = old_sensors(index);
    old_c_matrix(sensor_location,sensor_location) = 1;
end
old_og_matlab = lyap(total_a',old_c_matrix'*total_m*old_c_matrix,[],total_m');
old_og_matlab2 = lyap((total_m\eye(n+2,n+2))*total_a,old_c_matrix*old_c_matrix);


total_m_inv = total_m\eye(n+2,n+2);
a_minus_mu = (total_m_inv*total_a-mu*eye(n+2,n+2));
a_plus_mu = (total_m_inv*total_a+mu*eye(n+2,n+2));
a_plus_mu_inv = a_plus_mu\eye(n+2,n+2);

a_mu = a_minus_mu*a_plus_mu_inv;
c_mu_helper_left = (total_m_inv*total_a+mu*eye(n+2,n+2))\eye(n+2,n+2);
c_mu_helper_right = ((total_m_inv*total_a)'+mu*eye(n+2,n+2))\eye(n+2,n+2);

%diagonalize_a_mu
[a_mu_v,a_mu_d] = eig(a_mu);
a_mu_v_inv = a_mu_v\eye(n+2,n+2);

%new sensor index
new_sensor = 3;
new_c_matrix = zeros(n+2,n+2);
new_c_matrix(new_sensor,new_sensor) = 1;

%sylvester determinant theorem
sylvester_a = [];
sylvester_b = [];

%zeroth term
c_rank_1 = c_mu_helper_left*new_c_matrix;
c_rank_1_trans = new_c_matrix*c_mu_helper_right;
sylvester_a = c_mu_helper_left(:,new_sensor);
sylvester_b = c_mu_helper_right(new_sensor,:);


%for loop for the infinity series expansion
v_inv_c_rank_1 = a_mu_v_inv*c_rank_1;
c_rank_1_trans_v = c_rank_1_trans*(a_mu_v_inv)';

for term_index = 1:max_term
    
    sylvester_a_col = a_mu_v*(a_mu_d^term_index)*v_inv_c_rank_1;
    sylvester_a_col = sylvester_a_col(:,new_sensor);
    sylvester_a = [sylvester_a,sylvester_a_col]; %stack column by column
    
    
    sylvester_b_row = c_rank_1_trans_v*(a_mu_d^term_index)*a_mu_v';
    sylvester_b_row = sylvester_b_row(new_sensor,:);
    sylvester_b = [sylvester_b;sylvester_b_row]; %stack row by row
    
end

%now, we can compute the full gramian

new_c_total_matrix = old_c_matrix+new_c_matrix;
new_og_matlab2 = lyap((total_m\eye(n+2,n+2))*total_a,new_c_total_matrix*new_c_total_matrix);

diff = new_og_matlab2-(old_og_matlab2+(-2*mu)*sylvester_a*sylvester_a');



outputArg1 = norm(diff,'fro');


%check the difference of the determinant
[~,s_matlab,~] = svd(new_og_matlab2);
s_matlab = s_matlab(1:20,1:20);
log_det_matlab = sum(log(diag(s_matlab)));

[~,s_series,~] = svd(old_og_matlab2+(-2*mu)*sylvester_a*sylvester_a');
s_series = s_series(1:20,1:20);
log_det_series = sum(log(diag(s_series)));

outputArg2 = abs(log_det_matlab-log_det_series);

end

