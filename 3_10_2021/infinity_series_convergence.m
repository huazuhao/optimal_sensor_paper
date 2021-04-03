
clc;
clear;

fro_norm_diff = [];
objective_diff = [];

n = 100;
mu = -1;
max_terms = 100;

for max_term = 0:max_terms
    [outputArg1,outputArg2] = infinity_series_diff_function(n,mu,max_term);
    fro_norm_diff = [fro_norm_diff,outputArg1];
    objective_diff = [objective_diff,outputArg2];
end

f = figure;
plot(fro_norm_diff','linewidth',3);
%legend('full access to data','exact method','approximate method');
legend('frobenius norm diff between exact and infinity series for n = 100');
ylabel('frobenius norm')
xlabel('max term count')

figure;
plot(objective_diff','linewidth',3);
legend('difference of log det objective');
ylabel('log det objective')
xlabel('max term count')