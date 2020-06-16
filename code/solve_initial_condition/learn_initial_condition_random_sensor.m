clc;
clear;

rng(1);


%%

% Problem parameters
n     = 300;  % Number of interior spatial dofs
N     = 300;  % Number of temporal dofs
theta = 0.5;  % Trapezoidal rule, this parameter is used for discretization
gamma = 1e-1; % Robin coefficient, this is part of the boundry condition
alpha = 0.0000;    % Control penalty parameter, this parameter originates from the objective function
heat_constant = 1; %The constant in the heat equation. 

% Initial condition
x  = linspace(0, 1, n+2).';
t  = linspace(0, 1, N+2).';
U0 = 20*sin(6*pi*x)+cos(pi*x);
noise = normrnd(0,2,[n+2,1]);
U0 = U0+noise;
%plot(U0)
left_end = U0(1);
right_end = U0(n+2);
Z = rand(2,N+2); 
Z(1,:) = 20*sin(0.2*pi*t)+left_end;
Z(2,:) = 2*sin(0.3*pi*t)+right_end;

% Build finite element matrices
[M, A, B, R] = buildFEM(n, gamma);
A = heat_constant*A;
B = heat_constant*B;
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
random_selected_obs = position_index(1:10); %40 percent observation
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
% solve optimization problem with hessian times vector

state.is_state_computed = false;
state.is_adjoint_computed = false;

gtol = 1e-12;
U_init_recover = newton_method(gtol, state);


figure,
plot(x,U_init_recover,'r-','LineWidth',3), hold on
plot(x,U0,'k--','LineWidth',3), hold off
legend('Computed-with hessian times vector','True','Location','NorthEast')
xlabel('x','FontSize',20)
ylabel('U_0','FontSize',20)
set(gca,'FontSize',16)
%print('-depsc2','results-hess.eps')


%plot difference
figure
diff = U0-U_init_recover;
plot(diff);