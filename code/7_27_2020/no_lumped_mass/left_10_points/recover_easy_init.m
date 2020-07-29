clc;
clear;

rng(4);


%%
%values I am recording
norm_diff = [];
og_observation = [];




%%
%setting up the problem


% Problem parameters
n     = 100;  % Number of interior spatial dofs
N     = 100;  % Number of temporal dofs
theta = 0.5;  % Trapezoidal rule, this parameter is used for discretization
gamma = 1e-1; % Robin coefficient, this is part of the boundry condition
alpha = 0.0001; % Control penalty parameter, this parameter originates from the objective function
heat_constant = 1; %The constant in the heat equation. 
new_basis_rank = 10;
num_of_random_sensor = 5;

% Initial condition
x  = linspace(0, 1, n+2).';
t  = linspace(0, 1, N+2).';
U0 = 20*sin(8*pi*x)+cos(2*pi*x);
%U0 = 20*sin(4*pi*x)+cos(pi*x);
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


%compute the projection matrix

new_basis_rank = 10;
new_basis = zeros(n+2,new_basis_rank);
for index=1:new_basis_rank
    new_basis(:,index) = total_m(:,index);
end

new_basis = inv(total_m)*new_basis;
projection_matrix = inv(new_basis'*new_basis)*new_basis';

%first we are going to recover the initial condition based on algorithm
%selection
saved_matrix = matfile('approximation_method.mat');
saved_matrix = matfile('exact_method_10_choose_5.mat');
%saved_matrix = matfile('exact_method_10_choose_10.mat');
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
random_diff = U0_no_noise_projected-random_projected;
subspace_mass_matrix = total_m(1:new_basis_rank,1:new_basis_rank);
algorithm_diff = algorithm_diff'*subspace_mass_matrix*algorithm_diff
random_diff = random_diff'*subspace_mass_matrix*random_diff
