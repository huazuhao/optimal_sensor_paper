%%% BUILD FINITE ELEMENT MATRICES
function [M, A, B, R] = buildFEM(n, gamma)

% Define mesh size
h = 1/(n+1);

% Build mass matrix
% The first appearance of the mass matrix is on page 5
% Equation 1.3
% Since the basis function we are using is the hat function
% The realization of the mass matrix in 1d with the given basis function
% takes this form
o = h/6*ones(n+2, 1);
d = 4*h/6*ones(n+2, 1);
d(1)   = 2*h/6;
d(n+2) = 2*h/6;
M = spdiags([o d o], -1:1, n+2, n+2);

% Build stiffness matrix
% The first appearance of the mass matrix is on page 5
% equation 1.3
% This A matrix is usually called the stiffness matrix
o = -1/h*ones(n+2, 1);
d = 2/h*ones(n+2, 1);
d(1)   = 1/h + gamma;
d(n+2) = 1/h + gamma;
A = spdiags([o d o], -1:1, n+2, n+2);

% Build control operator
% This first appeared on page 5, equation 1.3
% Given the basis function, this is B
B = sparse(n+2, 2);
B(1,1)   = gamma;
B(n+2,2) = gamma;

% Build control penalty matrix
% This first sppeared on page 5, equation 1.3
% Given the basis function, this is R
R = speye(2,2);

return