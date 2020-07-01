clc;
clear;

rng(0);


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
alpha = 0.00001; % Control penalty parameter, this parameter originates from the objective function
heat_constant = 1; %The constant in the heat equation. 

% Initial condition
x  = linspace(0, 1, n+2).';
t  = linspace(0, 1, N+2).';
U0 = 20*sin(8*pi*x)+cos(2*pi*x);
%U0 = 20*sin(4*pi*x)+cos(pi*x);
noise = normrnd(0,2,[n+2,1]);
U0_no_noise = U0;
U0 = U0+noise;
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

best_recovery_norm = norm(U0_no_noise-U_init_recover,'fro');


%%
%here, we want to build a for loop

selected_obs = []; %this records what points are being selected for observation

for iteration_index = 1:10
    
    disp(iteration_index)
    
    %random access
    state.is_state_computed = false;
    state.is_adjoint_computed = false;
    position_index = randperm(n+2);
    random_selected_obs = position_index(1:iteration_index);
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
    
    
    
    random_norm = norm(U0_no_noise-U_init_recover_random_access,'fro');
    
    %i also want to automatically save images
    f = figure('visible','off');
    plot(U_init_recover_random_access,'r-','LineWidth',3);
    hold on;
    plot(U0,'k--','LineWidth',3)
    legend('recovered','truth');
    xlabel('x','FontSize',20)
    ylabel('U_0','FontSize',20)
    print(f,'-djpeg',sprintf('random_figure%d.jpg',iteration_index));
    delete(f);
    
    
    %greedy algorithm for maximizing og
    c_matrix = zeros(n+2,n+2);
    for index = 1:length(selected_obs)
        non_zero_index = selected_obs(index);
        c_matrix(non_zero_index, non_zero_index) = 1;
    end
    next_potential_placement = ones(1,n+2);
    for index = 1:length(selected_obs)
        already_selected = selected_obs(index);
        next_potential_placement(already_selected) = 0;
    end
    
    log_det_og_matlab = ones(1,n+2)*(-1e16);
    
    %parfor index = 1:n+2
    for index = 1:n+2    
        %construct a new c matrix
        if next_potential_placement(index) ~= 0
            new_c_matrix = c_matrix;
            new_c_matrix(index,index) = 1;
            og_matlab = lyap(transpose(total_a),new_c_matrix'*total_m*new_c_matrix,[],transpose(total_m));
            %og_matlab = lyap(transpose(total_a),new_c_matrix'*new_c_matrix,[],transpose(total_m));
            [~,s_matlab,~] = svd(og_matlab);
            log_det_matlab = sum(log(diag(s_matlab)));
            log_det_og_matlab(index) = log_det_matlab;
        end
    end
    [~,optimal_sensor_placement] = max(log_det_og_matlab);
    selected_obs = [selected_obs,optimal_sensor_placement];
    all_index = zeros(n+2,1);
    all_index(selected_obs) = 1;
    og_observation = [og_observation,all_index];
    disp(selected_obs)
    
    state.is_state_computed = false;
    state.is_adjoint_computed = false;
    state.obs_op= selected_obs;
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
    
    
    og_norm = norm(U0_no_noise-U_init_recover_og_access,'fro');
    
    diff_vector = [best_recovery_norm;random_norm;og_norm];
    diff_vector = [random_norm;og_norm];
    norm_diff = [norm_diff,diff_vector];
    
    %i also want to automatically save images
    f = figure('visible','off');
    plot(U_init_recover_og_access,'r-','LineWidth',3);
    hold on;
    plot(U0,'k--','LineWidth',3)
    legend('recovered','truth');
    xlabel('x','FontSize',20)
    ylabel('U_0','FontSize',20)
    print(f,'-djpeg',sprintf('og_figure%d.jpg',iteration_index));
    delete(f);
    
    
end


f = figure;
plot(norm_diff','linewidth',3);
%legend('full access to data','random access to data','og access to data');
legend('random access to data','og access to data');
ylabel('frobenius norm of difference between truth and recovered')
xlabel('sensor count')
print(f,'-djpeg',sprintf('result.jpg'));