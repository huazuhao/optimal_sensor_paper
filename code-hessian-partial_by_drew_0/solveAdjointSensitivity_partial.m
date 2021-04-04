function [P] = solveAdjointSensitivity_partial(S, state)

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

P = zeros(n+2, N+2);
% each column of P stores a solution to a time step of the adjoint sensitivity equation

% Recall Gmat = M - (1-theta)*dt*A
% Recall LU = (M+theta*dt*A)

Y = Lmat \ (theta*dt*Mobs*S(:,N+2));
% again, we start with the last time step for computing the equation. 
P(:,N+2) = Umat \ Y;
%the above two lines only solves the adjoint equation for the last time
%step

for i = N+1:-1:2

    Y = Lmat \ (Gmat*P(:,i+1) + dt*Mobs*S(:,i));
    P(:,i) = Umat \ Y;
    
end
% we have finished solving the adjoint sensitivity for the interior timing step. 

% solve for the first time
P(:,1) = Gmat*P(:,2)+dt*(1-theta)*Mobs*S(:,1);

return
