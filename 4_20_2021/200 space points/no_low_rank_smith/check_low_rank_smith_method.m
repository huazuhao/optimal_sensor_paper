clc;
clear;

rng(0);

%%
%what I mean by symmetric is that when I select the first sensor and the
%last sensor, the logdet of the observability gramian should be the same. 

%%
%setting up the problem


% Problem parameters
n     = 200;  % Number of interior spatial dofs
N     = 200;  % Number of temporal dofs
theta = 0.5;  % Trapezoidal rule, this parameter is used for discretization
gamma = 0.5; % Robin coefficient, this is part of the boundry condition
alpha = 0.00001; % Control penalty parameter, this parameter originates from the objective function
heat_constant = 1; %The constant in the heat equation. 
mu = -10; %mu, which should a number less than zero, from equation 3.3 of the paper titled "A modified low-rank smith method for large-scale lyapunov equations"
iteration_step_low_rank_smith = 100; %Number of iterations in the low rank smith method for computing the observability gramian



% Build finite element matrices
[M, A, B, R] = buildFEM(n, gamma);
A = heat_constant*A;
B = heat_constant*B;
total_a = full(-1*A);
total_m = full(M);
dt           = 1/(N+1);
[Lmat, Umat] = lu(M + theta*dt*A);
Gmat         = M - (1-theta)*dt*A;

%%
%low rank smith method
total_m_inv = total_m\eye(n+2,n+2);
total_m_inv_a = total_m_inv*total_a; %for the low rank smith method


selected_obs = [1,n+2,floor(n/2+1)]; %this records what points are being selected for observation
%old sensor locations
old_c_matrix = zeros(n+2,n+2);
for index = 1:length(selected_obs)
    sensor_location = selected_obs(index);
    old_c_matrix(sensor_location,sensor_location) = 1;
end

[og_pre_eig_vec,og_pre_eig_val] = low_rank_smith(mu,iteration_step_low_rank_smith,old_c_matrix,total_m_inv_a,n); %og_pre_eig_vec*og_pre_eig_val*eig_vec' = exact_og

og_low_rank_smith = og_pre_eig_vec*og_pre_eig_val*og_pre_eig_vec';

%%
%exact method

og_matlab = lyap(total_a',total_m'*old_c_matrix'*old_c_matrix*total_m,[],total_m');

%%
diff = norm(og_matlab-og_low_rank_smith,'fro')