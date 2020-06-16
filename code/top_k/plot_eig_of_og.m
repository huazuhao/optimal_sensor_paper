clc;
clear;

rng(0);



% Problem parameters
n     = 300;  % Number of interior spatial dofs
N     = 300;  % Number of temporal dofs
gamma = 0.5; % Robin coefficient, this is part of the boundary condition

[M, A, B, R] = buildFEM(n, gamma);

total_a = full(-1*A);
total_m = full(M);

selected_obs = randi([0,n+2],[1,60]);

new_c_matrix = zeros(n+2,n+2);
for index_2 = 1:length(selected_obs)
    non_zero_index = selected_obs(index_2);
    new_c_matrix(non_zero_index, non_zero_index) = 1;
end

og_matlab = lyap(transpose(total_a),new_c_matrix'*total_m*new_c_matrix,[],transpose(total_m));

[v,l] = eig(og_matlab);
[~,ind] = sort(diag(l),'descend');
l = l(ind,ind);
v = v(:,ind);

[~,s_matlab,~] = svd(og_matlab);

figure
plot(log(diag(s_matlab)),'LineWidth',5);
%plot(log(diag(l)),'LineWidth',5);
xlim([1,n+2])
xlabel('index')
ylabel('log of eigenvalue')
legend('eigenvalues of OG')

figure
plot_cut_off = 20;
selected_eig = diag(s_matlab);
%selected_eig = diag(l);
selected_eig = selected_eig(1:plot_cut_off);
plot(log(selected_eig),'LineWidth',5);
xlim([1,plot_cut_off])
xlabel('index')
ylabel('log of eigenvalue')
legend('eigenvalues of OG')

