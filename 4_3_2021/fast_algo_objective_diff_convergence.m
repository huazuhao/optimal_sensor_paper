
clc;
clear;

fro_norm_diff = [];
objective_diff = [];

n = 100;
mu = -1;
max_terms = 30;

for max_term = 0:max_terms
    
    disp(max_term)
    
    [outputArg] = fast_algo_objective_diff(max_term);
    objective_diff = [objective_diff,outputArg];
end


figure;
plot(objective_diff','linewidth',3);
legend('difference of log det objective');
ylabel('log det objective')
xlabel('max term count')