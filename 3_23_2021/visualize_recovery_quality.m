clc;
clear;

saved_matrix = matfile('recovery_quality_with_naive_implementation.mat');
naive_og = saved_matrix.norm_diff;

saved_matrix = matfile('recovery_quality_with_random.mat');
random_og = saved_matrix.norm_diff;

saved_matrix = matfile('recovery_quality_with_infinite_series_terms_0.mat');
series_0 = saved_matrix.norm_diff;
saved_matrix = matfile('recovery_quality_with_infinite_series_terms_10.mat');
series_10 = saved_matrix.norm_diff;
saved_matrix = matfile('recovery_quality_with_infinite_series_terms_100.mat');
series_100 = saved_matrix.norm_diff;

total_norm_diff = [series_0;series_10];
total_norm_diff = [random_og;naive_og;series_0;series_10];


f = figure;
plot(total_norm_diff','linewidth',3);
%legend('full access to data','random access to data','og access to data');
legend('random','exact','0th','10th');
ylabel('frobenius norm of difference between truth and recovered')
xlabel('sensor count')
