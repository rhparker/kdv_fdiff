% generate initial data
% takes too long, so we will load it from a file
% N = 1000;                 % number of finite diff gridpts
% iteration = 150;          % number of iterations of continuation code
% [x, uc] = solveKdV_fdiff(N, iterations);

% instead we load initial data from file
% this contains 150 iterations with step 0.25 in c
%  x:  grid of x values, 0 to 50, 1000 grid points
%  uc: continuation data; each column of uc is of the form
%       [u; c], where u is the solution and c is the speed
load uc_data

% specify which index you want in uc
index = 50
osc(x, uc, index)