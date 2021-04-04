function HTV = hessian_times_vector_partial(vector, U_init, state)

% Solve state equation
if (~state.is_state_computed)
  U = solveState(U_init, state);
else
  U = state.U;
end


% Solve adjoint equation
if (~state.is_adjoint_computed)
  L = solveAdjoint_partial(U, state);
else
  L = state.L;
end

% Solve sensitivity equation
S = solveSensitivity_partial(vector, state);

% Solve adjoint sensitivity equation
P = solveAdjointSensitivity_partial(S, state);

%return hessian times vector
HTV = -P(:,1);

return
