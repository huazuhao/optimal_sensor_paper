%%% SOLVE STATE EQUATION GIVEN CONTROL Z
function [U] = solveState(U_init, state)

B     = state.B;
Lmat  = state.Lmat;
Umat  = state.Umat;
Gmat  = state.Gmat;
N     = state.N;
n     = state.n;
dt    = state.dt;
theta = state.theta;
Z = rand(2,N+2); 
Z(1,:) = state.Z_left;
Z(2,:) = state.Z_right;

% Recall LU = (M+theta*dt*A)
% Recall Gmat = M - (1-theta)*dt*A
% we are trying to solve  LU x= Gmat+(dt*B*(theta*Z+(1-theta)*Z))
% This equation is the third to last equation on page 5 of Drew's original
% paper

U = zeros(n+2, N+2);
U(:,1) = U_init;
for i = 1:N+1
  Y = Lmat \ (Gmat*U(:,i) + dt*B*(theta*Z(:,i+1)+(1-theta)*Z(:,i)));
  U(:,i+1) = Umat \ Y;
end

% at the end of this program,
% we are returning a matrix u that is the solution of the PDE. 

return