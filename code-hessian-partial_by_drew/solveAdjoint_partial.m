function [L] = solveAdjoint(U, state)

M                  = state.M;
Mobs               = state.Mobs;
Lmat               = state.Lmat;
Umat               = state.Umat;
Gmat               = state.Gmat;
w                  = state.w;
N                  = state.N;
n                  = state.n;
dt                 = state.dt;
theta              = state.theta;

diff = Mobs*(U-w);

L = zeros(n+2, N+2);
% each column of L stores a solution to a time step of the adjoint equation

% Recall the adjoint equation
% (M+theta*dt*A)*lambda(:,i) = (-theta*dt*M*diff(:,i))
% Recall LU = (M+theta*dt*A)

Y = Lmat \ (-theta*dt*diff(:,N+2));
% again, we start with the last time step for computing the adjoint
% equation. 

L(:,N+2) = Umat \ Y;
%the above two lines only solves the adjoint equation for the last time
%step
%this is the same as drew's implementation

for i = N+1:-1:2
  Y = Lmat \ (Gmat*L(:,i+1) - dt*diff(:,i));
  L(:,i) = Umat \ Y;
end
% we have finished solving the adjoint for the interior timing step. 

L(:,1) = Gmat*L(:,2)-dt*(1-theta)*diff(:,1);

return

