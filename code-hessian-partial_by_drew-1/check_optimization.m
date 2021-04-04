clc;
clear;

rng(1);


%%

% Problem parameters
n     = 100;  % Number of interior spatial dofs
N     = 100;  % Number of temporal dofs
theta = 0.5;  % Trapezoidal rule, this parameter is used for discretization
gamma = 1e-1; % Robin coefficient, this is part of the boundry condition
alpha = 0;    % Control penalty parameter, this parameter originates from the objective function
heat_constant = 1; %The constant in the heat equation. 

% Initial condition
x  = linspace(0, 1, n+2).';
t  = linspace(0, 1, N+2).';
U0 = 20*sin(4*pi*x)+cos(pi*x);
left_end = U0(1);
right_end = U0(n+2);
Z = rand(2,N+2); 
Z(1,:) = 20*sin(0.2*pi*t)+left_end;
Z(2,:) = 2*sin(0.3*pi*t)+right_end;

% Build finite element matrices
[M, A, B, R] = buildFem(n, gamma);
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


%random access for spatial sensor
position_index = randperm(n+2);
random_selected_obs = position_index(1:5); %40 percent observation
state.obs_op= random_selected_obs; %obs_op means observation operator

Mobs = M;
for i = 1:n+2
  if i~=state.obs_op
    Mobs(i,:) = 0;
    Mobs(:,i) = 0;
  end
end

state.Mobs = Mobs;

%% 
% check gradient and Hessian computaitons
vector = rand(n+2,1);
random_initial_guess = rand(n+2,1);
[f, grad] = reducedObj_partial(random_initial_guess, state);
HTV = hessian_times_vector_partial(vector, random_initial_guess, state);
gv = grad.'*vector;

fprintf('  Check Gradient\n')
fprintf('  stepsize       dir deriv      finite diff    error\n')
eta = 1;
for i = 1:12
  xv   = random_initial_guess+eta*vector;
  fnew = reducedObj_partial(xv, state);
  fd   = (fnew - f)/eta;
  err  = abs(fd - gv);
  fprintf('  % 8.6e  % 8.6e  % 8.6e  % 8.6e\n',eta,gv,fd,err);
  eta  = eta*0.1;
end
fprintf('\n')

fprintf('  Check Hessian Times A Vector\n')
fprintf('  stepsize       hess vec       finite diff    error\n')
eta = 1;
for i = 1:12
    xv = random_initial_guess+eta*vector;
    [~,gxv] = reducedObj_partial(xv, state);
    fd = (gxv - grad)/eta;
    err = norm(fd - HTV);
    fprintf('  % 8.6e  % 8.6e  % 8.6e  % 8.6e\n',eta,norm(HTV),norm(fd),err);
    eta = eta*0.1;
end
fprintf('\n')

state.is_state_computed = false;
state.is_adjoint_computed = false;


%%
% now, we are going to solve the optimization problem with the old first
% derivative problem

% U_guess   = zeros(n+2,1);
% 
% options = optimoptions(@fminunc, ...
%                      'Algorithm', 'trust-region', ...
%                      'Display', 'iter', ...
%                      'SpecifyObjectiveGradient', true, ...
%                      'MaxIterations',100,...
%                      'FunctionToleranc',1e-12,...
%                      'OptimalityTolerance',1e-12);
%                  
% [U_init_recover,fmin_info,exitflag,output] = ...
%    fminunc(@(U_init)reducedObj_partial(U_init, state), U_guess, options);
% 
% 
% figure,
% plot(x,U_init_recover,'r-','LineWidth',3), hold on
% plot(x,U0,'k--','LineWidth',3), hold off
% legend('Computed-with only grad','True','Location','NorthEast')
% xlabel('x','FontSize',20)
% ylabel('U_0','FontSize',20)
% set(gca,'FontSize',16)
% print('-depsc2','results-nohess.eps')



%%
% solve optimization problem with hessian times vector

state.is_state_computed = false;
state.is_adjoint_computed = false;

%U_guess   = U_init_recover;
U_guess   = zeros(n+2,1);
[~, grad] = reducedObj_partial(U_guess, state);
grad      = -grad;
U_init_recover = pcg(@(v)hessian_times_vector_partial(v,U_guess,state),grad,1e-18,n+2);
U_init_recover = U_init_recover+U_guess;

% H = zeros(n+2,n+2);
% for i = 1:n+2
%   v = zeros(n+2,1); v(i) = 1;
%   H(:,i) = hessian_times_vector_partial(v,U_guess,state);
% end
% U_init_recover = H \ grad;
% U_init_recover = U_init_recover+U_guess;

f = reducedObj_partial(U_init_recover,state);
fprintf('  Final Objective Value: % 8.6e\n',f)


figure,
plot(x,U_init_recover,'r-','LineWidth',3), hold on
plot(x,U0,'k--','LineWidth',3), hold off
legend('Computed-with hessian times vector','True','Location','NorthEast')
xlabel('x','FontSize',20)
ylabel('U_0','FontSize',20)
set(gca,'FontSize',16)
print('-depsc2','results-hess.eps')
