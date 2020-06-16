%saved_matrix = matfile('approximation_10_method.mat');
saved_matrix = matfile('approximation_method.mat');
low_rank = saved_matrix.all_index;


%saved_matrix = matfile('exact_10_method.mat');
saved_matrix = matfile('exact_method.mat');
full_rank = saved_matrix.all_index;

figure,
plot(full_rank,'LineWidth',5);
hold on;
plot(low_rank,'*','LineWidth',5);

xlabel('point index')
ylabel('whether observe')
legend('low rank from current method','low rank from matlab')

%use top 10 eigenvalue
%we use 70 as the integration rank
%we are looking at the first 40 sensors

%use top 8 eigenvalue
%we use 30 as the integration rank
%we are looking at the first 30 sensors