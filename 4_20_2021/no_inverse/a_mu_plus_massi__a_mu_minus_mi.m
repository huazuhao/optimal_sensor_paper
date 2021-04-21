clc;
clear;

rng(0);

%%
%what I mean by symmetric is that when I select the first sensor and the
%last sensor, the logdet of the observability gramian should be the same. 

%%
%values I am recording
norm_diff = [];
%%
%setting up the problem


% Problem parameters
n     = 10;  % Number of interior spatial dofs
N     = 10;  % Number of temporal dofs
theta = 0.5;  % Trapezoidal rule, this parameter is used for discretization
gamma = 0.5; % Robin coefficient, this is part of the boundry condition
alpha = 0.00001; % Control penalty parameter, this parameter originates from the objective function
heat_constant = 1; %The constant in the heat equation. 
mu = -10; %mu, which should a number less than zero, from equation 3.3 of the paper titled "A modified low-rank smith method for large-scale lyapunov equations"
max_term = 20; %max term of the infinity series expansion of the gramian. 
iteration_step_low_rank_smith = 50; %Number of iterations in the low rank smith method for computing the observability gramian
projection_rank_cut_off = n+2; %Set the dimensionality of the subspace spanned by the observability gramian
regularization = 1e-6; %regularization for inverting the observability gramian of the previous iteration

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


%common terms of the observability gramian infinity series expansion
a_plus_mu_mass = total_a+mu*total_m;
a_minus_mu_mass = total_a-mu*total_m;

[a_plus_mu_mass_eigvec,a_plus_mu_mass_eigval] = eig(a_plus_mu_mass); %A = VDV^{-1}
[a_minus_mu_mass_eigvec,a_minus_mu_mass_eigval] = eig(a_minus_mu_mass);


%%
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

L = (total_m_inv*total_a+mu*eye(n+2,n+2))\eye(n+2,n+2);
R = ((total_m_inv*total_a)'+mu*eye(n+2,n+2))\eye(n+2,n+2);

[l_eigvec,l_eigval] = eig(L);
[r_eigvec,r_eigval] = eig(R);