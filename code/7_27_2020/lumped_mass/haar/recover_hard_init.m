clc;
clear;

rng(4);



%%
%setting up the problem

% Problem parameters
n     = 100;  % Number of interior spatial dofs
N     = 100;  % Number of temporal dofs
theta = 0.5;  % Trapezoidal rule, this parameter is used for discretization
gamma = 1e-1; % Robin coefficient, this is part of the boundry condition
alpha = 0.0001; % Control penalty parameter, this parameter originates from the objective function
heat_constant = 1; %The constant in the heat equation. 
new_basis_rank = 6;
num_of_random_sensor = 10;

% Initial condition
x  = linspace(0, 1, n+2).';
t  = linspace(0, 1, N+2).';
U0 = zeros(n+2,1);
U0(1:30) = 20;
U0(30:50) = -20;
U0(80:n+2) = -20;
noise = normrnd(0,2,[n+2,1]);
U0_no_noise = U0;
%U0 = U0+noise;
left_end = U0(1);
right_end = U0(n+2);
Z = rand(2,N+2); 
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


%%

%lumped mass matrix
lumped_m = zeros(n+2,n+2);
for index=1:n+2
    diag_sum = sum(total_m(index,:)); %row sum
    lumped_m(index,index) = diag_sum;
end

%compute a projection matrix 
saved_matrix = matfile('exact_haar_index_6.mat');
haar_index = saved_matrix.haar_index;


selected_j = [];
selected_k = [];
new_basis_vectors = [];

for index=1:new_basis_rank
    
    haar_j = haar_index(1,index);
    haar_k = haar_index(2,index);
    
    a_vector = haar_hat_function( haar_j,haar_k,n );
    
    new_basis_vectors = [new_basis_vectors,a_vector];
    
end

new_basis = inv(lumped_m)*new_basis_vectors;
projection_matrix = inv(new_basis'*new_basis)*new_basis';



%%

%first we are going to recover the initial condition based on algorithm
%selection
saved_matrix = matfile('approximation_method.mat');
saved_matrix = matfile('exact_method_6_choose_10.mat');
%saved_matrix = matfile('exact_method_5_choose_10.mat');
algorithm_sensors = saved_matrix.all_index;
selected_obs = [];

for index = 1:n+2
    if algorithm_sensors(index)==1
        selected_obs = [selected_obs,index];
    end
end
state.U = zeros(size(U));
state.L = zeros(size(U));
state.is_state_computed = false;
state.is_adjoint_computed = false;
state.obs_op = selected_obs;
Mobs = M;
for i = 1:n+2
  if i~=state.obs_op
    Mobs(i,:) = 0;
    Mobs(:,i) = 0;
  end
end
state.Mobs = Mobs;
gtol = 1e-12;
U_init_recover_algorithm_access = newton_method(gtol, state);


algorithm_projected = projection_matrix*U_init_recover_algorithm_access;


%random observation
state.is_state_computed = false;
state.is_adjoint_computed = false;
position_index = randperm(n+2);
random_selected_obs = position_index(1:num_of_random_sensor);
state.obs_op= random_selected_obs; %obs_op means observation operator
Mobs = M;
for i = 1:n+2
  if i~=state.obs_op
    Mobs(i,:) = 0;
    Mobs(:,i) = 0;
  end
end
state.Mobs = Mobs;
gtol = 1e-12;
U_init_recover_random_access = newton_method(gtol, state);

random_projected = projection_matrix*U_init_recover_random_access;

figure,
plot(U_init_recover_algorithm_access,'g','LineWidth',3);
hold on;
plot(U_init_recover_random_access,'k','LineWidth',3);
hold on;
plot(U0_no_noise,'LineWidth',3);
xlabel('position')
ylabel('temperature profile')
scatter(random_selected_obs,U_init_recover_random_access(random_selected_obs),80,'o','filled');
scatter(selected_obs,U_init_recover_algorithm_access(selected_obs),80,'d','filled');
legend('algorithm','random','ground truth','random sensor','algorithm sensor')


U0_no_noise_projected = projection_matrix*U0_no_noise;
% figure,
% plot(low_rank_projected,'g','LineWidth',3);
% hold on;
% plot(random_projected,'k','LineWidth',3);
% hold on;
% plot(U0_no_noise_projected,'LineWidth',3);
% legend('algorithm','random','ground truth')
% xlabel('position')
% ylabel('temperature profile in the selected subspace')


algorithm_diff = U0_no_noise_projected-algorithm_projected;
algorithm_diff = algorithm_diff'*algorithm_diff
random_diff = U0_no_noise_projected-random_projected;
random_diff = random_diff'*random_diff


