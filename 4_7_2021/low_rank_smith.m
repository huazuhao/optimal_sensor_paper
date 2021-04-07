function [eig_vec,eig_val] = low_rank_smith(mu,iteration_step,current_c,total_a,n)

    %implement the low rank smith paper
    u = [mu]; %the shift coefficients
    kshifts = 1;
    ksteps = iteration_step;

    Bm = current_c;
    Res = [];
    for j = 1:kshifts
        mu = u(j);
        rho = sqrt(-1*2*mu);
        Bm = ((total_a+mu*eye(n+2))\Bm);
        Res = [Res rho*Bm];
        for k = 1:ksteps
            Bm = (total_a-mu*eye(n+2))*((total_a+mu*eye(n+2))\Bm);
            Res = [Res rho*Bm];
        end
        Bm = (total_a-mu*eye(n+2))*Bm;
    end
    [V,S,Q] = svd(Res,0);

    eig_vec = V;
    eig_val = S*S';

end

