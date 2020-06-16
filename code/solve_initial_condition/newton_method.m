
function U_init_recover = newton_method(gtol, state)

%now, we will implement newton-cg with armigo line-search
continue_index_k = true;
uk = zeros(state.n+2,1); %initial guess
k_counter = 0;

while continue_index_k == true
    
    %step 2, compute gradient
    [objective_at_uk,grad] = reducedObj_partial(uk, state);
    
    %step 3, check stopping criteria
    if norm(grad)<gtol  %2-norm
        continue_index_k = false;
    end
    
    %step 5, implement the CG loop
    eta = 0.5;
    sk = 0;
    pki = -grad;
    rk0 = -grad;
    rki = rk0;
    continue_index_i = true;
    i_counter = -1;
    while continue_index_i == true
        
        i_counter = i_counter+1;
        
        %1. 
        if norm(rki,2)<eta*norm(rk0,2)
            continue_index_i = false;
        end
        %2.
        qki = hessian_times_vector_partial(pki,uk,state);
        %3.
        if pki'*qki<0
            continue_index_i = false;
        end
        %4.
        gammaki = norm(rki)^2/(pki'*qki);
        %5.
        sk = sk + gammaki*pki;
        %6
        rkiplus1 = rki-gammaki*qki;
        %7
        betaki = norm(rkiplus1)^2/norm(rki)^2;
        rki = rkiplus1;
        %8
        pki = rki + betaki*pki;
        
     
    end
    %step 5.3
    if i_counter == 0 
        sk = -grad;
    end
    
    %step 6
    %armigo line-search
    alpha_k = 1;
    new_objective_value = reducedObj_partial(uk+alpha_k*sk, state);  
    %only one output so the gradient is not required to compute
    while new_objective_value>objective_at_uk + 10e-4*alpha_k*sk'*grad
        alpha_k = alpha_k/2;
        new_objective_value = reducedObj_partial(uk+alpha_k*sk, state);
    end
   
    %step 7
    uk = uk+alpha_k*sk;
    
    
    k_counter = k_counter+1;
    disp(k_counter)
    if k_counter>=50
        continue_index_k = false;
    end
    
end

U_init_recover = uk;