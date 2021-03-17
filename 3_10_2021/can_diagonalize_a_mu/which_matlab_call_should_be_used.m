clc;
clear;


% Problem parameters
n     = 10;  % Number of interior spatial dofs
N     = 10;  % Number of temporal dofs
theta = 0.5;  % Trapezoidal rule, this parameter is used for discretization
gamma = 0.5; % Robin coefficient, this is part of the boundry condition
alpha = 0.0001; % Control penalty parameter, this parameter originates from the objective function
heat_constant = 1; %The constant in the heat equation. 


[M, A, B, R] = buildFEM(n, gamma);
A = heat_constant*A;
B = heat_constant*B;
total_a = full(-1*A);
total_m = full(M);


%given the old sensor and compute the observability gramian
old_sensors = [1,5,9];
old_c_matrix = zeros(n+2,n+2);
for index = 1:length(old_sensors)
    sensor_location = old_sensors(index);
    old_c_matrix(sensor_location,sensor_location) = 1;
end
old_og_matlab = lyap(total_a',old_c_matrix'*total_m*old_c_matrix,[],total_m');
old_og_matlab2 = lyap((total_m\eye(n+2,n+2))*total_a,old_c_matrix*old_c_matrix);


norm(old_og_matlab-old_og_matlab2,'fro')