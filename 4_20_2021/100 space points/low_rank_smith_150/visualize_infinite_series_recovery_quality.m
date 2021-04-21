clc;
clear;

saved_matrix = matfile('recovery_quality_with_infinite_series_terms_0.mat');
series_0 = saved_matrix.norm_diff;
saved_matrix = matfile('recovery_quality_with_infinite_series_terms_10.mat');
series_10 = saved_matrix.norm_diff;
saved_matrix = matfile('recovery_quality_with_infinite_series_terms_20.mat');
series_20 = saved_matrix.norm_diff;
%saved_matrix = matfile('recovery_quality_with_infinite_series_terms_100.mat');
%series_100 = saved_matrix.norm_diff;


f = figure;
x_tick = linspace(1,13,13);


x_tick = linspace(4,13,10);
plot(x_tick,series_0,'linewidth',3);
hold on;
plot(x_tick,series_10,'linewidth',3);
hold on;
plot(x_tick,series_20,'linewidth',3);
hold on;


legend('0th','10th','20th');
ylabel('frobenius norm of difference between truth and recovered')
xlabel('sensor count (low rank smith with 150 iterations) ')
xlim([1,13])

