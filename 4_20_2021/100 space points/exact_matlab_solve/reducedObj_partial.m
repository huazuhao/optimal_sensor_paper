function [val, grad] = reducedObj(U_init, state)
% Solve state equation
U = solveState(U_init, state);

% Store computed state
state.U = U;
state.is_state_computed = true;

% Compute objective value
M                  = state.M;
Mobs               = state.Mobs;
w                  = state.w;
N                  = state.N;
n                  = state.n;
alpha              = state.alpha;
theta              = state.theta;
dt                 = state.dt;
Gmat               = state.Gmat;

diff = U-w; % the difference between the current state and the final objective
Mdiff = Mobs*diff;
MU_init = M*U_init;

% compute the objective function
val = (1-theta)*(diff(:,1).'*Mdiff(:,1)); %beginning of time
for i =1:N
    val = val + diff(:,i+1).'*Mdiff(:,i+1);
end
val = val + theta*(diff(:,N+2).'*Mdiff(:,N+2)); %end of time
val = val*dt;
val = val + alpha*U_init.'*MU_init; % since we are guessing initial condition now
val = val*0.5; 

if nargout > 1 % if grad is also requested to compute
  % Solve adjoint equation
  L = solveAdjoint_partial(U, state);

  % Store computed adjoint
  state.L = L;
  state.is_adjoint_computed = true;

  % Build gradient
  grad = alpha*MU_init-L(:,1);

end

return
