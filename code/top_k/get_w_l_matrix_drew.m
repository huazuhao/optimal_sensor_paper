function [w,l] = get_w_l_matrix_drew(M,A)


%%
A = A*-1;


%%
%now, I am going to perform those four steps described by Drew

%step 1
[u,k] = eig(full(M)); %i could use cholesky and then only look at the diagonal entries. 

%step 2
u_hat = u*inv(k^(1/2));
a_hat = u_hat'*A*u_hat;

%step 3
[v,l] = eig(a_hat);
%sort v and l so that the largest eigenvalue is at the end of l
[~,ind] = sort(diag(l),'descend');
%[~,ind] = sort(diag(l));

l = l(ind,ind);
v = v(:,ind);

%step 4
w = u_hat*v;
