clc;
clear;

saved_matrix = matfile('recovery_quality_with_naive_implementation.mat');
naive_og = saved_matrix.norm_diff;


saved_matrix = matfile('recovery_quality_average_of_random.mat');
experiment_data = saved_matrix.experiment_data;
random_og = mean(experiment_data,1);

random_mean = mean(experiment_data,1);
random_std = std(experiment_data,0,1);

saved_matrix = matfile('recovery_quality_with_infinite_series_terms_0.mat');
series_0 = saved_matrix.norm_diff;
saved_matrix = matfile('recovery_quality_with_infinite_series_terms_10.mat');
series_10 = saved_matrix.norm_diff;
saved_matrix = matfile('recovery_quality_with_infinite_series_terms_30.mat');
series_30 = saved_matrix.norm_diff;
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
plot(x_tick,series_30,'linewidth',3);
hold on;
% plot(x_tick,series_50,'linewidth',3);
% hold on;

%legend('full access to data','random access to data','og access to data');
legend('random','exact','0th','10th','30th');
ylabel('frobenius norm of difference between truth and recovered')
xlabel('sensor count')
xlim([1,13])

