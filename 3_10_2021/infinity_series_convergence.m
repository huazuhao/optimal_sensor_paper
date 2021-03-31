
clc;
clear;

fro_norm_diff = [];
n = 100;
mu = -1;
max_terms = 20;

for max_term = 0:max_terms
    [outputArg] = infinity_series_diff_function(n,mu,max_term);
    fro_norm_diff = [fro_norm_diff,outputArg];
end

f = figure;
plot(fro_norm_diff','linewidth',3);
%legend('full access to data','exact method','approximate method');
legend('frobenius norm diff between exact and infinity series for n = 100');
ylabel('frobenius norm')
xlabel('max term count')