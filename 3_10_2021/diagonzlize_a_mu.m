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


[M, A, B, R] = buildFEM(n, gamma);
A = heat_constant*A;
B = heat_constant*B;
total_a = full(-1*A);
total_m = full(M);

total_m_inv = total_m\eye(n+2,n+2);
a_minus_mu = (total_m_inv*total_a-mu*eye(n+2,n+2));
a_plus_mu = (total_m_inv*total_a+mu*eye(n+2,n+2));
a_plus_mu_inv = a_plus_mu\eye(n+2,n+2);

a_mu = a_minus_mu*a_plus_mu_inv;
b_mu_helper_left = (total_m_inv*total_a+mu*eye(n+2,n+2))\eye(n+2,n+2);
b_mu_helper_right = ((total_m_inv*total_a)'+mu*eye(n+2,n+2))\eye(n+2,n+2);


diff = a_mu-a_mu';
diff_2 = b_mu_helper_left-b_mu_helper_right; 

%now, we diagonalize a_mu
[a_mu_v,a_mu_d] = eig(a_mu);

diff_3 = a_mu-a_mu_v*a_mu_d*(a_mu_v\eye(n+2,n+2));

diff_4 = a_mu*a_mu-a_mu_v*(a_mu_d*a_mu_d)*(a_mu_v\eye(n+2,n+2));
