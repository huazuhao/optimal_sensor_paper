function [W] = solveAdjoint(vector, state)

M                  = state.M;
Lmat               = state.Lmat;
Umat               = state.Umat;
Gmat               = state.Gmat;
w                  = state.w;
N                  = state.N;
n                  = state.n;
dt                 = state.dt;
theta              = state.theta;
alpha              = state.alpha;


W = zeros(n+2, N+2);
% each column of W stores a solution to a time step of the sensitivity equation

% Recall Gmat = M - (1-theta)*dt*A
% Recall LU = (M+theta*dt*A)


W(:,1) = -vector;

for i = 2:N+2
    Y = Lmat \ (Gmat * W(:,i-1));
    W(:,i) = Umat \ Y;
end

return

