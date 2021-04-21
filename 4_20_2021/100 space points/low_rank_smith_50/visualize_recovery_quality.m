clc;
clear;

% Initial condition
n = 100;
N = 100;
x  = linspace(0, 1, n+2).';
t  = linspace(0, 1, N+2).';
U0 = 20*sin(8*pi*x)+cos(2*pi*x);
initial_condition_norm = norm(U0,'fro');

saved_matrix = matfile('recovery_quality_with_naive_implementation.mat');
naive_og = saved_matrix.norm_diff;
naive_og = naive_og/initial_condition_norm;

saved_matrix = matfile('recovery_quality_average_of_random.mat');
experiment_data = saved_matrix.experiment_data;
experiment_data = experiment_data/initial_condition_norm;
random_og = mean(experiment_data,1);

random_mean = mean(experiment_data,1);
random_std = std(experiment_data,0,1);

saved_matrix = matfile('recovery_quality_with_infinite_series_terms_0.mat');
series_0 = saved_matrix.norm_diff;
series_0 = series_0/initial_condition_norm;
saved_matrix = matfile('recovery_quality_with_infinite_series_terms_10.mat');
series_10 = saved_matrix.norm_diff;
series_10 = series_10/initial_condition_norm;
saved_matrix = matfile('recovery_quality_with_infinite_series_terms_20.mat');
series_20 = saved_matrix.norm_diff;
series_20 = series_20/initial_condition_norm;
% saved_matrix = matfile('recovery_quality_with_infinite_series_terms_50.mat');
% series_50 = saved_matrix.norm_diff;


f = figure;
x_tick = linspace(1,13,13);

errorbar(x_tick,random_mean,random_std,'LineWidth',5);
hold on;

plot(x_tick,naive_og,'linewidth',3);
hold on;

x_tick = linspace(4,13,10);
plot(x_tick,series_0,'linewidth',3);
hold on;
plot(x_tick,series_10,'linewidth',3);
hold on;
plot(x_tick,series_20,'linewidth',3);
hold on;
% plot(x_tick,series_50,'linewidth',3);
% hold on;

%legend('full access to data','random access to data','og access to data');
legend('random','exact','0th','10th','20th');
ylabel('fro norm of dif between truth and recovered/fro norm of truth')
xlabel('sensor count (low rank smith 50 iterations)')
xlim([1,13])

