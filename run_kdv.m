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

%% load data and make configuration
load uc_data;
load uc_full;

% method to use
config.method   = 'fdiff';
config.BC       = 'Neumann';

% specify which index you want in uc
index       = 35;
index_full  = 48;

% input data to use
uin     = uc(:,index);
u_full  = uc_full(:,index_full);

% % ramp up to 5000 grid points
% % don't think we actually need to do this
% [xin, uin, c] = solveKdV_fdiff_newton_interp(x, u, config, 5000);

%% make plot of oscillations
% plot oscillations and get scaled version of u
c = uin(end);
[y, uscaled, start, decay, freq] = osc_plot_c(x, uin);
% scaled version with same domain as uin
% can plot this alongside uin
uscaled2 = [ zeros(start-1,1); uscaled ];

%% find the minima and maxima of this
udata = uin(1:end-1);
uscaled = udata.*exp(decay*x);

% compute finite difference matrices
N = length(x);
L = x(end) - x(1);             
h = L/(N-1); 
D = D_fdiff(N, h, config.BC);
y = D*uscaled./uscaled - decay;

% threshold for a zero of this function
epsilon = 0.025;
zeros = find(abs(y) < epsilon);

%% use output of half-wave to make full wave
[x2, u2] = full_wave(x, uin);

% Newton solver using Neumann BCs
un_out = solveKdV_fdiff_newton(x2, u2, config);

% Newton solver using periodic BCs
p_config.method = 'fdiff'; 
p_config.BC     = 'periodic';
% need to remove last point for periodic BCs
% since it is identified with first point
xper = x2(1:end-1);
uper = [u2(1:end-2); c];
up_out = solveKdV_fdiff_newton(xper, uper, p_config);

plot(x2, un_out(1:end-1), xper, up_out(1:end-1))

%% eigenvalues of linearization
 