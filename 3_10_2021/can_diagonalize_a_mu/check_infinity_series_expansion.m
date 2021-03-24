clc;
clear;


% Problem parameters
n     = 10;  % Number of interior spatial dofs
N     = 10;  % Number of temporal dofs
theta = 0.5;  % Trapezoidal rule, this parameter is used for discretization
gamma = 0.5; % Robin coefficient, this is part of the boundry condition
alpha = 0.0001; % Control penalty parameter, this parameter originates from the objective function
heat_constant = 1; %The constant in the heat equation. 
mu = -1;
max_term = 3; %max term of the infinity series expansion of the gramian. 


[M, A, B, R] = buildFEM(n, gamma);
A = heat_constant*A;
B = heat_constant*B;
total_a = full(-1*A);
total_m = full(M);


%given the old sensor and compute the observability gramian
old_sensors = [];
old_c_matrix = zeros(n+2,n+2);
for index = 1:length(old_sensors)
    sensor_location = old_sensors(index);
    old_c_matrix(sensor_location,sensor_location) = 1;
end
old_og_matlab = lyap(total_a',total_m'*(old_c_matrix'*old_c_matrix)*total_m,[],total_m');
old_og_matlab2 = lyap(total_m\total_a,old_c_matrix*old_c_matrix);


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
new_sensor = n+2;
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
new_og_infinity_series = old_og_matlab + (-2*mu)*sylvester_a*sylvester_b

new_c_total_matrix = old_c_matrix+new_c_matrix;
new_og_matlab = lyap(total_a',total_m'*(new_c_total_matrix'*new_c_total_matrix)*total_m,[],total_m');
new_og_matlab2 = lyap((total_m\eye(n+2,n+2))*total_a,new_c_total_matrix*new_c_total_matrix);

new_og_matlab

diff = new_og_matlab-new_og_infinity_series;
diff = new_og_matlab2-(old_og_matlab2+(-2*mu)*sylvester_a*sylvester_b);


