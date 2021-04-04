
function [outputArg] = fast_algo_objective_diff(max_term)


rng(0);

%%
%what I mean by symmetric is that when I select the first sensor and the
%last sensor, the logdet of the observability gramian should be the same.
%%
%setting up the problem


% Problem parameters
n     = 100;  % Number of interior spatial dofs
N     = 100;  % Number of temporal dofs
theta = 0.5;  % Trapezoidal rule, this parameter is used for discretization
gamma = 0.5; % Robin coefficient, this is part of the boundry condition
alpha = 0.00001; % Control penalty parameter, this parameter originates from the objective function
heat_constant = 1; %The constant in the heat equation. 
mu = -10; %mu, which should a number less than zero, from equation 3.3 of the paper titled "A modified low-rank smith method for large-scale lyapunov equations"
iteration_step_low_rank_smith = 50; %Number of iterations in the low rank smith method for computing the observability gramian
projection_rank_cut_off = 20; %Set the dimensionality of the subspace spanned by the observability gramian

% Initial condition
x  = linspace(0, 1, n+2).';
t  = linspace(0, 1, N+2).';
U0 = 20*sin(8*pi*x)+cos(2*pi*x);
%U0 = 20*sin(4*pi*x)+cos(pi*x);
left_end = U0(1);
right_end = U0(n+2);
Z = rand(2,N+2); %Z contains boundary conditions
Z(1,:) = 20*sin(0.2*pi*t)+left_end;
Z(2,:) = 2*sin(0.3*pi*t)+right_end;


% Build finite element matrices
[M, A, B, R] = buildFEM(n, gamma);
A = heat_constant*A;
B = heat_constant*B;
total_a = full(-1*A);
total_m = full(M);
dt           = 1/(N+1);
[Lmat, Umat] = lu(M + theta*dt*A);
Gmat         = M - (1-theta)*dt*A;

% Set state memory
state.n     = n;
state.N     = N;
state.theta = theta;
state.gamma = gamma;
state.alpha = alpha;
state.M     = M;
state.A     = A;
state.B     = B;
state.R     = R;
state.dt    = dt;
state.Lmat  = Lmat;
state.Umat  = Umat;
state.Gmat  = Gmat;
state.u0    = U0;
state.Z_left  = Z(1,:);
state.Z_right = Z(2,:);

U = solveState(U0, state);
state.w = U;

state.U = zeros(size(U));
state.L = zeros(size(U));
state.is_state_computed = false;
state.is_adjoint_computed = false;


%common terms of the observability gramian infinity series expansion
total_m_inv = total_m\eye(n+2,n+2);
total_m_inv_a = total_m_inv*total_a; %for the low rank smith method
a_minus_mu = (total_m_inv*total_a-mu*eye(n+2,n+2));
a_plus_mu = (total_m_inv*total_a+mu*eye(n+2,n+2));
a_plus_mu_inv = a_plus_mu\eye(n+2,n+2);
a_mu = a_minus_mu*a_plus_mu_inv;
%diagonalize_a_mu
[a_mu_v,a_mu_d] = eig(a_mu);
a_mu_v_inv = a_mu_v\eye(n+2,n+2);

c_mu_helper_left = (total_m_inv*total_a+mu*eye(n+2,n+2))\eye(n+2,n+2);
c_mu_helper_right = ((total_m_inv*total_a)'+mu*eye(n+2,n+2))\eye(n+2,n+2);
%c_mu_helper_right = c_mu_helper_left';

lyapunov_inf_series.a_mu_v = a_mu_v;
lyapunov_inf_series.a_mu_d = a_mu_d;
lyapunov_inf_series.a_mu_v_inv = a_mu_v_inv;
lyapunov_inf_series.c_mu_helper_left = c_mu_helper_left;
lyapunov_inf_series.c_mu_helper_right = c_mu_helper_right;

%%
%the first step we want to do is to establish the best boundary recovery we
%can do with full access to data
%recover initial condition with full access to data
% state.Mobs = M;
% gtol = 1e-12;
% U_init_recover = newton_method(gtol, state);
% best_recovery_norm = norm(U0-U_init_recover,'fro');


%%
%we now begin to select sensors. 

selected_obs = [1,n+2,floor(n/2+1)]; %this records what points are being selected for observation


%old sensor locations
old_c_matrix = zeros(n+2,n+2);
for index = 1:length(selected_obs)
    sensor_location = selected_obs(index);
    old_c_matrix(sensor_location,sensor_location) = 1;
end


%old_og_matlab = lyap(total_a',total_m'*old_c_matrix'*old_c_matrix*total_m,[],total_m');
%[og_pre_eig_vec,og_pre_eig_val] = eig(old_og_matlab);
[og_pre_eig_vec,og_pre_eig_val] = low_rank_smith(mu,iteration_step_low_rank_smith,old_c_matrix,total_m_inv_a,n); %og_pre_eig_vec*og_pre_eig_val*eig_vec' = exact_og

index = 10; %position of the new sensor
new_c_matrix = zeros(n+2,n+2);
new_c_matrix(index,index) = 1;


%sylvester determinant theorem
sylvester_a = [];
sylvester_b = [];
%zeroth term
sylvester_a = lyapunov_inf_series.c_mu_helper_left(:,index);
sylvester_b = lyapunov_inf_series.c_mu_helper_right(index,:);
%common terms in the infinite series
c_rank_1 = lyapunov_inf_series.c_mu_helper_left*new_c_matrix;
c_rank_1_trans = new_c_matrix*lyapunov_inf_series.c_mu_helper_right;
v_inv_c_rank_1 = lyapunov_inf_series.a_mu_v_inv*c_rank_1;
c_rank_1_trans_v = c_rank_1_trans*(lyapunov_inf_series.a_mu_v_inv)';

for term_index = 1:max_term
    sylvester_a_col = lyapunov_inf_series.a_mu_v*(lyapunov_inf_series.a_mu_d^term_index)*v_inv_c_rank_1;
    sylvester_a_col = sylvester_a_col(:,index);
    sylvester_a = [sylvester_a,sylvester_a_col]; %stack column by column


    sylvester_b_row = c_rank_1_trans_v*(lyapunov_inf_series.a_mu_d^term_index)*(lyapunov_inf_series.a_mu_v)';
    sylvester_b_row = sylvester_b_row(index,:);
    sylvester_b = [sylvester_b;sylvester_b_row]; %stack row by row

end

%now, compute the log det objective
sigma_pre_top = og_pre_eig_val(1:projection_rank_cut_off,1:projection_rank_cut_off);
sigma_pre_top_inv = sigma_pre_top\eye(projection_rank_cut_off,projection_rank_cut_off);

u_pre = og_pre_eig_vec(:,1:projection_rank_cut_off);
objective_det_sigma = sum(log(diag(sigma_pre_top)));
objective_det_rest = eye(max_term+1,max_term+1)-2*mu*sylvester_a'*u_pre*sigma_pre_top_inv*u_pre'*sylvester_a;

[~,s_next,~] = svd(objective_det_rest);
log_det_obj = sum(log(diag(s_next)));



%exact from matlab
total_c_matrix = old_c_matrix+new_c_matrix;
og_matlab = lyap(total_a',total_m'*total_c_matrix'*total_c_matrix*total_m,[],total_m');
[og_next_eig_vec,og_next_eig_val] = eig(og_matlab);
[og_next_eig_val,sorted_ind] = sort(diag(og_next_eig_val),'descend');
og_next_eig_vec = og_next_eig_vec(:,sorted_ind);
log_det_obj_matlab = sum(log(og_next_eig_val(1:projection_rank_cut_off)));

outputArg = abs(objective_det_sigma+log_det_obj-log_det_obj_matlab);

%without projection
% [~,s_series,~] = svd(og_pre_eig_vec*og_pre_eig_val*og_pre_eig_vec'+(-2*mu)*sylvester_a*sylvester_a');
% log_det_series = sum(log(diag(s_series)));
% 
% [~,s_matlab,~] = svd(og_matlab);
% log_det_matlab = sum(log(diag(s_matlab)));
% outputArg = abs(log_det_series-log_det_matlab);



%diff = og_matlab-(og_pre_eig_vec*og_pre_eig_val*og_pre_eig_vec'+(-2*mu)*sylvester_a*sylvester_a');
%outputArg = norm(diff,'fro');
