clc;
clear;

rng(0);

%%
%values I am recording
norm_diff = [];
%%
%setting up the problem


% Problem parameters
n     = 100;  % Number of interior spatial dofs
N     = 100;  % Number of temporal dofs
theta = 0.5;  % Trapezoidal rule, this parameter is used for discretization
gamma = 0.5; % Robin coefficient, this is part of the boundry condition
alpha = 0.00001; % Control penalty parameter, this parameter originates from the objective function
heat_constant = 1; %The constant in the heat equation. 

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

%%
%the first step we want to do is to establish the best boundary recovery we
%can do with full access to data
%recover initial condition with full access to data
state.Mobs = M;
gtol = 1e-12;
U_init_recover = newton_method(gtol, state);

best_recovery_norm = norm(U0-U_init_recover,'fro');


%%
%here, we want to build a for loop for selecting sensors

selected_obs = []; %this records what points are being selected for observation

for iteration_index = 1:10
    
    disp(iteration_index)
    
    
    %random access
    state.is_state_computed = false;
    state.is_adjoint_computed = false;
    position_index = randperm(n+2);
    random_selected_obs = position_index(1:iteration_index);
    state.obs_op= random_selected_obs; %obs_op means observation operator

    
    
    state.is_state_computed = false;
    state.is_adjoint_computed = false;
    Mobs = M;
    for i = 1:n+2
      if i~=state.obs_op
        Mobs(i,:) = 0;
        Mobs(:,i) = 0;
      end
    end
    state.Mobs = Mobs;
    gtol = 1e-12;
    U_init_recover_og_access = newton_method(gtol, state);
    
    series_norm = norm(U0-U_init_recover_og_access,'fro');
    
    diff_vector = [best_recovery_norm;series_norm];
    diff_vector = [series_norm];
    norm_diff = [norm_diff,diff_vector];
    
 
end

all_index = zeros(1,n+2);
all_index(selected_obs) = 1;
figure;
plot(all_index);

f = figure;
plot(norm_diff','linewidth',3);
%legend('full access to data','random access to data','og access to data');
legend('exact og access to data');
ylabel('frobenius norm of difference between truth and recovered')
xlabel('sensor count')

save('recovery_quality_with_random.mat','norm_diff')